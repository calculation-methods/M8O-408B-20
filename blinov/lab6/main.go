package main

import (
	"fmt"
	"image/color"
	"math"
)

func solution(x, t float64) float64 {
	return math.Cos(x) * math.Sin(2*t)
}

//# boundary conditions
//def phi_0(t):
//return math.sin(2*t)
//
//def phi_1(t):
//return -math.sin(2*t)
//
//# initial conditions
//def psi_0(x):
//return 0
//
//def psi_1(x):
//return 2 * math.cos(x)
//
//def solution(x, t):
//return math.cos(x) * math.sin(2*t)

func phi0(x float64) float64 {
	return math.Sin(2 * x)
}

func phiL(t float64) float64 {
	return -math.Sin(2 * t)
}

func psi0(x float64) float64 {
	return 0
}

func psiL(x float64) float64 {
	return 2 * math.Cos(x)
}

func analyticalSolution(
	xRange, tRange Interval,
	h, sigma float64,
) [][]float64 {
	tau := math.Sqrt(sigma * math.Pow(h, 2))
	x := xRange.Range(h)
	t := tRange.Range(tau)

	res := make([][]float64, len(t))
	for idt := range t {
		res[idt] = make([]float64, len(x))
		for idx, xVal := range x {
			res[idt][idx] = solution(xVal, t[idt])
		}
	}

	return res
}

func explicitFiniteDifferenceMethod(xRange, tRange Interval, h, sigma, a float64) [][]float64 {
	tau := math.Sqrt(sigma * math.Pow(h, 2)) // len of cell by t
	x := xRange.Range(h)
	t := tRange.Range(tau)

	res := make([][]float64, len(t))
	for i := range res {
		res[i] = make([]float64, len(x))
	}

	// Row 0 -> use initial condition 0
	for colId, xVal := range x {
		res[0][colId] = psi0(xVal)
	}

	// Row 1 -> use approximation
	for colId, xVal := range x {
		res[1][colId] = psi0(xVal) + tau*psiL(xVal)
	}

	for rowId := 2; rowId < len(t); rowId++ {
		// Col 0 -> use boundary condition 0
		res[rowId][0] = phi0(t[rowId])
		// Cols 1..n-1 -> use explicit schema
		for colId := 1; colId < len(x)-1; colId++ {
			res[rowId][colId] = sigma*(res[rowId-1][colId+1]-2*res[rowId-1][colId]+res[rowId-1][colId-1]) +
				(2-3*math.Pow(tau, 2))*res[rowId-1][colId] - res[rowId-2][colId]
		}
		// Col n -> use boundary condition 1
		res[rowId][len(x)-1] = phiL(t[rowId])
	}
	return res
}

func tridiagonalSolve(A [][]float64, b []float64) ([]float64, error) {
	n := len(A)
	if n == 0 || len(b) != n {
		return nil, fmt.Errorf("invalid dimensions for A or b")
	}

	v := make([]float64, n)
	u := make([]float64, n)
	v[0] = A[0][1] / -A[0][0]
	u[0] = b[0] / A[0][0]
	for i := 1; i < n-1; i++ {
		v[i] = A[i][i+1] / (-A[i][i] - A[i][i-1]*v[i-1])
		u[i] = (A[i][i-1]*u[i-1] - b[i]) / (-A[i][i] - A[i][i-1]*v[i-1])
	}
	v[n-1] = 0
	u[n-1] = (A[n-1][n-2]*u[n-2] - b[n-1]) / (-A[n-1][n-1] - A[n-1][n-2]*v[n-2])

	x := make([]float64, n)
	x[n-1] = u[n-1]
	for i := n - 1; i > 0; i-- {
		x[i-1] = v[i-1]*x[i] + u[i-1]
	}
	return x, nil
}

func fillTridiagonalMatrix(A [][]float64, sigma, tau float64) {
	A[0][0] = -(1 + 2*sigma + 3*math.Pow(tau, 2))
	A[0][1] = sigma
	for i := 1; i < len(A)-1; i++ {
		A[i][i-1] = sigma
		A[i][i] = -(1 + 2*sigma + 3*math.Pow(tau, 2))
		A[i][i+1] = sigma
	}
	A[len(A)-1][len(A)-2] = sigma
	A[len(A)-1][len(A)-1] = -(1 + 2*sigma + 3*math.Pow(tau, 2))
}

func implicitFiniteDifferenceMethod(
	xRange, tRange Interval,
	h, sigma, a float64,
) [][]float64 {
	tau := math.Sqrt(sigma * math.Pow(h, 2)) // len of cell by t
	x := xRange.Range(h)
	t := tRange.Range(tau)

	res := make([][]float64, len(t))
	for i := range res {
		res[i] = make([]float64, len(x))
	}

	for colId, xVal := range x {
		res[0][colId] = psi0(xVal)
		res[1][colId] = psi0(xVal) + tau*psiL(xVal)
	}

	for rowId := 2; rowId < len(t); rowId++ {
		n := len(x) - 2
		A := make([][]float64, n)
		for i := range A {
			A[i] = make([]float64, n)
		}

		fillTridiagonalMatrix(A, sigma, tau)

		b := make([]float64, n)
		for i := range b {
			b[i] = -2*res[rowId-1][i+1] + res[rowId-2][i+1]
		}
		b[0] -= sigma * phi0(t[rowId])
		b[n-1] -= sigma * phiL(t[rowId])

		solution, err := tridiagonalSolve(A, b)
		if err != nil {
			panic(fmt.Sprintf("error solving tridiagonal system: %v", err))
		}

		res[rowId][0] = phi0(t[rowId])
		res[rowId][len(x)-1] = phiL(t[rowId])
		copy(res[rowId][1:len(x)-1], solution)
	}

	return res
}

func main() {

	a := 1.

	xRange := Interval{0, math.Pi}
	tRange := Interval{0, 5}

	h := 0.01
	sigma := 1.

	analyticalSolutionResult := analyticalSolution(xRange, tRange, h, sigma)

	// Explicit
	fmt.Println("Explicit")
	explicitSolutionResult := explicitFiniteDifferenceMethod(xRange, tRange, h, sigma, a)
	maxAbsErr, _ := maxAbsError2D(analyticalSolutionResult, explicitSolutionResult)
	meanAbsErr, _ := meanAbsError(analyticalSolutionResult, explicitSolutionResult)
	fmt.Printf("max abs error = %.12f\n", maxAbsErr)
	fmt.Printf("mean abs error = %.12f\n", meanAbsErr)

	// Implicit
	fmt.Println("\nImplicit")
	implicitSolutionResult := implicitFiniteDifferenceMethod(xRange, tRange, h, sigma, a)
	maxAbsErr, _ = maxAbsError2D(analyticalSolutionResult, implicitSolutionResult)
	meanAbsErr, _ = meanAbsError(analyticalSolutionResult, implicitSolutionResult)
	fmt.Printf("max abs error = %.12f\n", maxAbsErr)
	fmt.Printf("mean abs error = %.12f\n", meanAbsErr)

	plotResults(map[string]struct {
		arr   [][]float64
		color color.Color
	}{
		"analytical": {analyticalSolutionResult, color.RGBA{255, 0, 0, 255}},
		"explicit":   {explicitSolutionResult, color.RGBA{0, 255, 0, 255}},
		"implicit":   {implicitSolutionResult, color.RGBA{0, 0, 255, 255}},
	}, 0.5, xRange, tRange, h, sigma, a)
	//
	plotErrorsFromTime(map[string]struct {
		arr   [][]float64
		color color.Color
	}{
		"analytical": {analyticalSolutionResult, color.RGBA{255, 0, 0, 255}},
		"explicit":   {explicitSolutionResult, color.RGBA{0, 255, 0, 255}},
		"implicit":   {implicitSolutionResult, color.RGBA{0, 0, 255, 255}},
	}, tRange, h, sigma, "error_plot1.png")

	plotErrorsFromTime(map[string]struct {
		arr   [][]float64
		color color.Color
	}{
		"analytical": {analyticalSolutionResult, color.RGBA{255, 0, 0, 255}},
		"explicit":   {explicitSolutionResult, color.RGBA{0, 255, 0, 255}},
	}, tRange, h, sigma, "error_plot2.png")
}
