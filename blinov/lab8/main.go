package main

import (
	"fmt"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
	"math"
)

func solution(ap, x, y, t float64) float64 {
	return math.Cos(2*x) * math.Cosh(y) * math.Exp(-3*ap*t)
}

func analytFunc(ap float64, X, Y [][]float64, t float64) [][]float64 {
	rows := len(Y)
	cols := len(X[0])
	result := make([][]float64, rows)
	for i := 0; i < rows; i++ {
		result[i] = make([]float64, cols)
		for j := 0; j < cols; j++ {
			result[i][j] = math.Exp(-3*ap*t) * math.Cos(2*X[i][j]) * math.Cosh(Y[i][j]) // Adjust the function as needed
		}
	}
	return result
}

func u0yt(ap, y, t float64) float64 {
	return math.Cosh(y) * math.Exp(-3*ap*t)
}

func upi4yt(ap, y, t float64) float64 {
	return 0
}

func ux0t(ap, x, t float64) float64 {
	return math.Cos(2*x) * math.Exp(-3*ap*t)
}

func uxln2t(ap, x, t float64) float64 {
	return 1.25 * math.Cos(2*x) * math.Exp(-3*ap*t)
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

func runThrough(a, b, c, d []float64, s int) []float64 {
	P := make([]float64, s+1)
	Q := make([]float64, s+1)

	P[0] = -c[0] / b[0]
	Q[0] = d[0] / b[0]

	k := s - 1
	for i := 1; i < s; i++ {
		denominator := b[i] + a[i]*P[i-1]
		P[i] = -c[i] / denominator
		Q[i] = (d[i] - a[i]*Q[i-1]) / denominator
	}
	P[k] = 0
	Q[k] = (d[k] - a[k]*Q[k-1]) / (b[k] + a[k]*P[k-1])

	x := make([]float64, s)
	x[k] = Q[k]

	for i := s - 2; i >= 0; i-- {
		x[i] = P[i]*x[i+1] + Q[i]
	}

	return x
}

// variableDirections implements the method of variable directions for solving PDEs.
func variableDirections(ap float64, x, y []float64, hx, hy float64, K int, tau, t float64) [][][]float64 {
	U := make([][][]float64, K)
	for k := range U {
		U[k] = make([][]float64, len(x))
		for i := range U[k] {
			U[k][i] = make([]float64, len(y))
		}
	}

	sigmaA := (ap * tau) / (2 * float64(hx*hx))

	// Set initial border conditions
	for k := 0; k < K; k++ {
		for i := range x {
			U[k][i][0] = ux0t(ap, x[i], t)
			U[k][i][len(y)-1] = uxln2t(ap, x[i], t)
		}
		for j := range y {
			U[k][0][j] = u0yt(ap, y[j], t)
			U[k][len(x)-1][j] = upi4yt(ap, y[j], t)
		}
	}

	// Main computation loop
	for k := 1; k < K; k++ {
		Utemp := make([][]float64, len(x))
		for i := range Utemp {
			Utemp[i] = make([]float64, len(y))
		}

		t += tau / 2
		// First half-step
		for i := 1; i < len(x)-1; i++ {
			a := make([]float64, len(y))
			b := make([]float64, len(y))
			c := make([]float64, len(y))
			d := make([]float64, len(y))

			for j := 1; j < len(y)-1; j++ {
				a[j] = sigmaA
				b[j] = -1 - 2*sigmaA
				c[j] = sigmaA
				d[j] = -sigmaA*(U[k-1][i][j+1]-2*U[k-1][i][j]+U[k-1][i][j-1]) + U[k-1][i][j]
			}

			// Boundary conditions
			b[0] = 1
			c[0] = 0
			d[0] = ux0t(ap, x[i], t-tau/2)

			a[len(y)-1] = 0
			b[len(y)-1] = 1
			d[len(y)-1] = uxln2t(ap, x[i], t-tau/2)

			uNew := runThrough(a, b, c, d, len(y))
			Utemp[i] = uNew
		}

		// Update boundary values
		for j := range y {
			Utemp[0][j] = u0yt(ap, y[j], t)
			Utemp[len(x)-1][j] = upi4yt(ap, y[j], t)
		}

		// Second half-step
		t += tau / 2
		for j := 1; j < len(y)-1; j++ {
			a := make([]float64, len(x))
			b := make([]float64, len(x))
			c := make([]float64, len(x))
			d := make([]float64, len(x))

			for i := 1; i < len(x)-1; i++ {
				a[i] = sigmaA
				b[i] = -1 - 2*sigmaA
				c[i] = sigmaA
				d[i] = sigmaA*(Utemp[i+1][j]-2*Utemp[i][j]+Utemp[i-1][j]) + Utemp[i][j]
			}

			// Boundary conditions
			b[0] = 1
			c[0] = 0
			d[0] = u0yt(ap, y[j], t)

			a[len(x)-1] = 0
			b[len(x)-1] = 1
			d[len(x)-1] = upi4yt(ap, y[j], t)

			uNew := runThrough(a, b, c, d, len(x))
			for i := range uNew {
				U[k][i][j] = uNew[i]
			}
		}

		// Update boundary values
		for i := range x {
			U[k][i][0] = ux0t(ap, x[i], t)
			U[k][i][len(y)-1] = uxln2t(ap, x[i], t)
		}

	}

	U = transpose3D(U)
	// Transpose U if needed
	return U
}

func fractionalStep(ap float64, x, y []float64, hx, hy float64, K int, tau, t float64) [][][]float64 {
	U := make([][][]float64, K)
	for k := range U {
		U[k] = make([][]float64, len(x))
		for i := range U[k] {
			U[k][i] = make([]float64, len(y))
		}
	}

	sigmaA := (ap * tau) / (hx * hx)

	// Initialize borders
	for k := 0; k < K; k++ {
		for i := range x {
			U[k][i][0] = ux0t(ap, x[i], t)
			U[k][i][len(y)-1] = uxln2t(ap, x[i], t)
		}
		for j := range y {
			U[k][0][j] = u0yt(ap, y[j], t)
			U[k][len(x)-1][j] = upi4yt(ap, y[j], t)
		}
	}

	// Main computation loop
	for k := 1; k < K; k++ {
		Utemp := make([][]float64, len(x))
		for i := range Utemp {
			Utemp[i] = make([]float64, len(y))
		}

		t += tau / 2
		// First part of the fractional step
		for i := 1; i < len(x)-1; i++ {
			a := make([]float64, len(y))
			b := make([]float64, len(y))
			c := make([]float64, len(y))
			d := make([]float64, len(y))

			for j := 1; j < len(y)-1; j++ {
				a[j] = sigmaA
				b[j] = -1 - 2*sigmaA
				c[j] = sigmaA
				d[j] = -Utemp[i][j]
			}

			// Boundary conditions
			b[0] = 1 - 0/hy
			c[0] = 0 / hy
			d[0] = ux0t(ap, x[i], t-tau/2)

			a[len(y)-1] = -0 / hy
			b[len(y)-1] = 1 + 0/hy
			d[len(y)-1] = uxln2t(ap, x[i], t-tau/2)

			uNew := runThrough(a, b, c, d, len(y))
			Utemp[i] = uNew
		}

		// Update boundary values
		for j := range y {
			Utemp[0][j] = u0yt(ap, y[j], t)
			Utemp[len(x)-1][j] = upi4yt(ap, y[j], t)
		}

		// Second part of the fractional step
		t += tau / 2
		for j := 1; j < len(y)-1; j++ {
			a := make([]float64, len(x))
			b := make([]float64, len(x))
			c := make([]float64, len(x))
			d := make([]float64, len(x))

			for i := 1; i < len(x)-1; i++ {
				a[i] = sigmaA
				b[i] = -1 - 2*sigmaA
				c[i] = sigmaA
				d[i] = -U[k-1][i][j]
			}

			// Boundary conditions
			b[0] = 1 - 0/hx
			c[0] = 0 / hx
			d[0] = u0yt(ap, y[j], t-tau/2)

			a[len(x)-1] = -0 / hx
			b[len(x)-1] = 1 + 0/hx
			d[len(x)-1] = upi4yt(ap, y[j], t-tau/2)

			uNew := runThrough(a, b, c, d, len(x))
			for i := range uNew {
				U[k][i][j] = uNew[i]
			}
		}

		// Update boundary values
		for i := range x {
			U[k][i][0] = ux0t(ap, x[i], t)
			U[k][i][len(y)-1] = uxln2t(ap, x[i], t)
		}
	}

	return transpose3D(U) // Assuming transpose3D function is implemented as previously described
}

func main() {
	// Inputs: Nx, Ny, K, time
	Nx := 30.
	Ny := 30.
	K := 300
	time := 10.

	hx := (math.Pi / 4) / float64(Nx)
	hy := math.Log(2) / float64(Ny)
	xrange := Interval{
		0, math.Pi / 4,
	}
	yrange := Interval{
		0, math.Log(2),
	}
	trange := Interval{
		0, time,
	}
	x := xrange.Range(hx)
	y := yrange.Range(hy)
	tau := time / float64(K)
	T := trange.Range(tau)
	t := 0.0

	for {
		fmt.Println("Выберите метод:\n1 - метод переменных направлений\n2 - метод дробных шагов\n0 - выход из программы")
		var method int
		fmt.Scanln(&method)
		if method == 0 {
			break
		}

		var ap float64
		fmt.Print("Введите параметр a -> ")
		fmt.Scanln(&ap)

		// variableDirections(ap float64, x, y []float64, hx, hy, tau float64, K int, t float64)

		var U [][][]float64
		if method == 1 {
			//fmt.Println(float64(ap), x, y, hx, hy, K, tau, t)
			U = variableDirections(float64(ap), x, y, hx, hy, K, tau, t)
		} else if method == 2 {
			U = fractionalStep(float64(ap), x, y, hx, hy, K, tau, t)
		}

		//
		var dt int
		fmt.Print("Введите момент времени -> ")
		fmt.Scanln(&dt)

		X := make([][]float64, len(y))
		Y := make([][]float64, len(x))
		for i := range X {
			X[i] = make([]float64, len(x))
			Y[i] = make([]float64, len(x))
			for j := range X[i] {
				X[i][j] = x[j]
				Y[i][j] = y[i]
			}
		}

		// Calculate U_analytic and errors
		U_analytic := make([][]float64, len(y))
		for i := range U_analytic {
			U_analytic[i] = make([]float64, len(x))
			for j := range U_analytic[i] {
				U_analytic[i][j] = solution(ap, x[j], y[i], T[dt])
			}
		}

		errorX := make([]float64, len(x))
		errorY := make([]float64, len(y))
		errorT := make([]float64, K)

		for i := range x {
			maxDiffY := 0.0
			for j := range y {
				diff := math.Abs(U_analytic[j][i] - U[dt][i][j])
				if diff > maxDiffY {
					maxDiffY = diff
				}
			}
			errorY[i] = maxDiffY
		}

		for j := range y {
			maxDiffX := 0.0
			for i := range x {
				diff := math.Abs(U_analytic[j][i] - U[dt][i][j])
				if diff > maxDiffX {
					maxDiffX = diff
				}
			}
			errorX[j] = maxDiffX
		}

		for k := 0; k < K; k++ {
			Uk := make([][]float64, len(y))
			for i := range Uk {
				Uk[i] = make([]float64, len(x))
				for j := range Uk[i] {
					if k < len(U) {
						Uk[i][j] = U[k][i][j]
					}
				}
			}
			errorT[k] = normalize(analytFunc(ap, X, Y, T[k]), Uk)
		}

		// Plotting errors (using gonum/plot)
		p := plot.New()

		p.Title.Text = "Error Graph"
		p.X.Label.Text = "x, y, t"
		p.Y.Label.Text = "error"

		// Add lines for each error type
		_ = plotutil.AddLinePoints(p,
			"Fixed x at a specific time", pointsFromSlice(x, errorY),
			"Fixed y at a specific time", pointsFromSlice(y, errorX),
			"All time intervals", pointsFromSlice(T, errorT))

		// Save the plot to a PNG file
		if err := p.Save(4*vg.Inch, 4*vg.Inch, "errors.png"); err != nil {
			panic(err)
		}

		fmt.Println("Error plot saved as 'errors.png'")
		break // Remove this if you want the loop to continue

	}

	return
}

func pointsFromSlice(x []float64, y []float64) plotter.XYs {
	pts := make(plotter.XYs, len(x))
	for i := range x {
		pts[i].X = x[i]
		pts[i].Y = y[i]
	}
	return pts
}
