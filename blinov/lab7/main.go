package main

import (
	"fmt"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
	"math"
)

const (
	ZEIDEL = iota
	LIEBMAB
	SIMPLE
)

func solution(x, y float64) float64 {
	return math.Exp(x) * math.Cos(y)
}

func u0y(y float64) float64 {
	return math.Cos(y)
}

func uly(y float64) float64 {
	return math.E * math.Cos(y)
}

func ux0(y float64) float64 {
	return 0
}

func uxl(x float64) float64 {
	return -math.Exp(x)
}

func normalize(current, previous [][]float64) float64 {
	maxDiff := 0.0
	for i := range current {
		for j := range current[i] {
			diff := math.Abs(current[i][j] - previous[i][j])
			if diff > maxDiff {
				maxDiff = diff
			}
		}
	}
	return maxDiff
}

func copyMatrix(dst, src [][]float64) {
	for i := range src {
		for j := range src[i] {
			dst[i][j] = src[i][j]
		}
	}
}

func liebman(x, y []float64, h [2]float64, eps float64) ([][]float64, int) {
	N := len(x)
	M := len(y)
	count := 0
	previous := make([][]float64, N)
	current := make([][]float64, N)
	for i := range previous {
		previous[i] = make([]float64, M)
		current[i] = make([]float64, M)
	}
	for j := range y {
		current[0][j] = u0y(y[j])
		current[N-1][j] = uly(y[j])
	}

	for j := range y {
		for i := 1; i < N-1; i++ {
			current[i][j] = current[i][0] + (current[i][M-1]-current[i][0])/(x[N-1]-x[0])*(x[i]-x[0])
		}
	}

	for normalize(current, previous) > eps {
		count++
		copyMatrix(previous, current)
		for i := 1; i < N-1; i++ {
			for j := 1; j < M-1; j++ {
				current[i][j] = (math.Pow(h[0], 2)*(previous[i-1][j]+previous[i+1][j]) +
					math.Pow(h[1], 2)*(previous[i][j-1]+previous[i][j+1])) / (2 * (math.Pow(h[0], 2) + math.Pow(h[1], 2)))
			}
		}
		for i := range x {
			current[i][0] = current[i][1] - h[1]*ux0(x[i])
			current[i][M-1] = current[i][M-2] + h[1]*uxl(x[i])
		}
	}

	U := make([][]float64, N)
	for i := 0; i < N; i++ {
		U[i] = make([]float64, M)
	}
	copyMatrix(U, current)
	return U, count
}

func relaxation(x, y []float64, h [2]float64, eps float64, w float64) ([][]float64, int) {
	N := len(x)
	M := len(y)
	count := 0

	previous := make([][]float64, N)
	current := make([][]float64, N)
	for i := range previous {
		previous[i] = make([]float64, M)
		current[i] = make([]float64, M)
	}

	for j := range y {
		current[0][j] = u0y(y[j])
		current[N-1][j] = uly(y[j])
	}

	for j := range y {
		for i := 1; i < N-1; i++ {
			current[i][j] = current[i][0] + (current[i][M-1]-current[i][0])/(x[N-1]-x[0])*(x[i]-x[0])
		}
	}

	for normalize(current, previous) > eps {
		count++
		copyMatrix(previous, current)
		for i := 1; i < N-1; i++ {
			for j := 1; j < M-1; j++ {
				current[i][j] = (math.Pow(h[0], 2)*(current[i-1][j]+previous[i+1][j]) +
					math.Pow(h[1], 2)*(current[i][j-1]+previous[i][j+1])) / (2 * (math.Pow(h[0], 2) + math.Pow(h[1], 2)))
				current[i][j] = current[i][j]*w + (1-w)*previous[i][j]
			}
		}

		for i := range x {
			current[i][0] = current[i][1] - h[1]*ux0(x[i])
			current[i][M-1] = current[i][M-2] + h[1]*uxl(x[i])
		}
	}

	U := make([][]float64, N)
	for i := 0; i < N; i++ {
		U[i] = make([]float64, M)
	}

	copyMatrix(U, current)
	return U, count
}

func Zeidel(x, y []float64, h [2]float64, eps float64) ([][]float64, int) {
	return relaxation(x, y, h, eps, 1) // вызываем метод релаксации с коэффициентом w=1
}

func main() {
	hx := 0.1
	hy := 0.1

	eps := 0.001
	//n := 10

	h := [2]float64{hx, hy}
	xrange := Interval{
		Start: 0,
		End:   1.05 + hx/2 - 1e-4,
	}

	yrange := &Interval{
		Start: 0,
		End:   math.Pi/2 + hy/2 - 1e-4,
	}

	x := xrange.Range(hx)
	y := yrange.Range(hy)

	for method := range []int{LIEBMAB, ZEIDEL, SIMPLE} {
		var U [][]float64
		var count int

		switch method {
		case LIEBMAB:
			U, count = liebman(x, y, h, eps)
		case ZEIDEL:
			U, count = Zeidel(x, y, h, eps)
		case SIMPLE:
			U, count = relaxation(x, y, h, eps, 1.8)
		}

		fmt.Println(U)

		fmt.Printf("Метод завершил работу с %d итерациями\n", count)
		// Здесь должен быть код для вывода U и сравнения с аналитическим решением

		// Расчет аналитического решения U_analytic и ошибок
		U_analytic := make([][]float64, len(x))
		for i := 0; i < len(x); i++ {
			U_analytic[i] = make([]float64, len(y))
		}

		errorX := make([]float64, len(x))
		errorY := make([]float64, len(y))
		for i, xv := range x {
			for j, yv := range y {
				if i == 0 {
					U_analytic[i] = make([]float64, len(y))
				}
				U_analytic[i][j] = solution(xv, yv)
				// Заполнение ошибок для каждого x и y
				diff := math.Abs(U_analytic[i][j] - U[i][j])
				if diff > errorX[i] {
					errorX[i] = diff
				}
				if diff > errorY[j] {
					errorY[j] = diff
				}
			}
		}

		// Создание нового графика.
		p := plot.New()

		// Настройка заголовка и осей графика.
		p.Title.Text = "График ошибок"
		p.X.Label.Text = "x, y"
		p.Y.Label.Text = "Ошибка"

		// Создание точек для графиков ошибок по x и y.
		pointsX := make(plotter.XYs, len(x))
		pointsY := make(plotter.XYs, len(y))
		for i := range x {
			pointsX[i].X = x[i]
			pointsX[i].Y = errorX[i]
		}
		for i := range y {
			pointsY[i].X = y[i]
			pointsY[i].Y = errorY[i]
		}

		// Добавление кривых ошибок на график.
		lineX, err := plotter.NewLine(pointsX)
		if err != nil {
			panic(err)
		}
		lineX.Color = plotutil.Color(1)

		lineY, err := plotter.NewLine(pointsY)
		if err != nil {
			panic(err)
		}
		lineY.Color = plotutil.Color(2)

		p.Add(lineX, lineY)
		p.Legend.Add("Ошибка при фиксированном x", lineX)
		p.Legend.Add("Ошибка при фиксированном y", lineY)

		// Сохранение графика в файл.
		if err := p.Save(8*vg.Inch, 4*vg.Inch, fmt.Sprintf("error_plot_%d.png", method)); err != nil {
			panic(err)
		}

		n := 3
		m := 3

		yStep := len(y) / (n * m)

		for i := 0; i < n; i++ {
			for j := 0; j < m; j++ {
				pY := j * yStep
				if pY >= len(y) {
					continue
				}

				p := plot.New()

				p.Title.Text = fmt.Sprintf("Сравнение решений при y = %.2f", y[pY])
				p.X.Label.Text = "x"
				p.Y.Label.Text = "u"

				// Добавление линий на график
				lineNum := make(plotter.XYs, len(x))
				lineAnalytic := make(plotter.XYs, len(x))
				for k, xVal := range x {
					lineNum[k].X = xVal
					lineNum[k].Y = U[k][pY]
					lineAnalytic[k].X = xVal
					lineAnalytic[k].Y = solution(xVal, y[pY])
				}

				numLine, _ := plotter.NewLine(lineNum)
				numLine.Color = plotutil.Color(0)

				analyticLine, _ := plotter.NewLine(lineAnalytic)
				analyticLine.Color = plotutil.Color(1)

				p.Add(numLine, analyticLine)
				p.Legend.Add("Численный метод", numLine)
				p.Legend.Add("Аналитическое решение", analyticLine)

				// Сохранение графика в файл
				fileName := fmt.Sprintf("%dplot_y_%.2f.png", method, y[pY])
				if err := p.Save(8*vg.Inch, 8*vg.Inch, fileName); err != nil {
					panic(err)
				}
			}
		}

		xStep := len(x) / (n * m)

		for i := 0; i < n; i++ {
			for j := 0; j < m; j++ {
				pX := i * xStep
				if pX >= len(x) {
					continue
				}

				p := plot.New()

				p.Title.Text = fmt.Sprintf("Сравнение решений при x = %.2f", x[pX])
				p.X.Label.Text = "y"
				p.Y.Label.Text = "u"

				// Добавление линий на график
				lineNum := make(plotter.XYs, len(y))
				lineAnalytic := make(plotter.XYs, len(y))
				for k, yVal := range y {
					lineNum[k].X = yVal
					lineNum[k].Y = U[pX][k]
					lineAnalytic[k].X = yVal
					lineAnalytic[k].Y = solution(x[pX], yVal)
				}

				numLine, _ := plotter.NewLine(lineNum)
				numLine.Color = plotutil.Color(0)

				analyticLine, _ := plotter.NewLine(lineAnalytic)
				analyticLine.Color = plotutil.Color(1)

				p.Add(numLine, analyticLine)
				p.Legend.Add("Численный метод", numLine)
				p.Legend.Add("Аналитическое решение", analyticLine)

				// Сохранение графика в файл
				fileName := fmt.Sprintf("%dplot_x_%.2f.png", method, x[pX])
				if err := p.Save(8*vg.Inch, 8*vg.Inch, fileName); err != nil {
					panic(err)
				}
			}
		}

	}
}
