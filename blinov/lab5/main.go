package main

import (
	"fmt"
	"image/color"
	"math"
)

func solution(x, t, a float64) float64 {
	return math.Exp(-a*t) * math.Cos(x)
}

func phi0(t, a float64) float64 {
	return math.Exp(-a * t)
}

func phiL(t, a float64) float64 {
	return -math.Exp(-a * t)
}

func psi(x float64) float64 {
	return math.Cos(x)
}

func analyticalSolution(xRange, tRange Interval, h, sigma, a float64) [][]float64 {
	tau := sigma * math.Pow(h, 2) / a

	x := xRange.Range(h)
	t := tRange.Range(tau)

	result := make([][]float64, len(t))
	for idt := range t {
		result[idt] = make([]float64, len(x))
		for idx := range x {
			result[idt][idx] = solution(x[idx], t[idt], a)
		}
	}

	return result
}

func explicitFiniteDifferenceMethod(xRange, tRange Interval, h, sigma, a float64) [][]float64 {
	tau := sigma * math.Pow(h, 2) / a
	x := xRange.Range(h)
	t := tRange.Range(tau)

	result := make([][]float64, len(t))
	for i := range result {
		result[i] = make([]float64, len(x))
	}

	for colId, xVal := range x {
		result[0][colId] = psi(xVal)
	}

	for rowId := 1; rowId < len(t); rowId++ {
		result[rowId][0] = phi0(t[rowId], a)

		for colId := 1; colId < len(x)-1; colId++ {
			result[rowId][colId] = sigma*result[rowId-1][colId-1] +
				(1-2*sigma)*result[rowId-1][colId] +
				sigma*result[rowId-1][colId+1]
		}

		result[rowId][len(x)-1] = phiL(t[rowId], a)
	}

	return result
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

func implicitFiniteDifferenceMethod(
	xRange, tRange Interval,
	h, sigma, a float64,
) ([][]float64, error) {
	tau := sigma * math.Pow(h, 2) / a
	x := xRange.Range(h)
	t := tRange.Range(tau)
	res := make([][]float64, len(t))
	for i := range res {
		res[i] = make([]float64, len(x))
	}

	for colId, xVal := range x {
		res[0][colId] = psi(xVal)
	}

	for rowId := 1; rowId < len(t); rowId++ {
		n := len(x) - 2
		A := make([][]float64, n)
		for i := range A {
			A[i] = make([]float64, n)
		}

		A[0][0] = -(1 + 2*sigma)
		A[0][1] = sigma
		for i := 1; i < n-1; i++ {
			A[i][i-1] = sigma
			A[i][i] = -(1 + 2*sigma)
			A[i][i+1] = sigma
		}
		A[n-1][n-2] = sigma
		A[n-1][n-1] = -(1 + 2*sigma)

		b := make([]float64, n)
		for i := range b {
			b[i] = -res[rowId-1][i+1]
		}

		b[0] -= sigma * phi0(t[rowId], a)
		b[n-1] -= sigma * phiL(t[rowId], a)

		slv, err := tridiagonalSolve(A, b)
		if err != nil {
			return nil, fmt.Errorf("error solving tridiagonal system: %v", err)
		}

		res[rowId][0] = phi0(t[rowId], a)
		res[rowId][len(x)-1] = phiL(t[rowId], a)
		copy(res[rowId][1:len(x)-1], slv)
	}

	return res, nil
}

func crankNicolsonMethod(
	xRange, tRange Interval,
	h, sigma, a, theta float64,
) ([][]float64, error) {
	tau := sigma * math.Pow(h, 2) / a
	x := xRange.Range(h)
	t := tRange.Range(tau)
	res := make([][]float64, len(t))
	for i := range res {
		res[i] = make([]float64, len(x))
	}

	// Row 0 -> use initial condition
	for colId, xVal := range x {
		res[0][colId] = psi(xVal)
	}

	for rowId := 1; rowId < len(t); rowId++ {
		n := len(x) - 2
		A := make([][]float64, n)
		for i := range A {
			A[i] = make([]float64, n)
		}

		// Create system of equations for implicit schema
		A[0][0] = -(1 + 2*sigma*theta)
		A[0][1] = sigma * theta
		for i := 1; i < n-1; i++ {
			A[i][i-1] = sigma * theta
			A[i][i] = -(1 + 2*sigma*theta)
			A[i][i+1] = sigma * theta
		}
		A[n-1][n-2] = sigma * theta
		A[n-1][n-1] = -(1 + 2*sigma*theta)

		// Vector b
		b := make([]float64, n)
		for i := range b {
			b[i] = -(res[rowId-1][i+1] + (1-theta)*sigma*(res[rowId-1][i]-2*res[rowId-1][i+1]+res[rowId-1][i+2]))
		}
		b[0] -= sigma * theta * phi0(t[rowId], a)
		b[n-1] -= sigma * theta * phiL(t[rowId], a)

		// Solve tridiagonal system
		solution, err := tridiagonalSolve(A, b)
		if err != nil {
			return nil, fmt.Errorf("error solving tridiagonal system: %v", err)
		}

		res[rowId][0] = phi0(t[rowId], a)
		res[rowId][len(x)-1] = phiL(t[rowId], a)
		copy(res[rowId][1:len(x)-1], solution)
	}

	return res, nil
}

func main() {

	a := 1.

	xRange := Interval{0, math.Pi}
	tRange := Interval{0, 5}

	h := 0.1
	sigma := 0.45

	analyticalSolutionResult := analyticalSolution(xRange, tRange, h, sigma, a)

	// Explicit
	fmt.Println("Explicit")
	explicitSolutionResult := explicitFiniteDifferenceMethod(xRange, tRange, h, sigma, a)
	maxAbsErr, _ := maxAbsError2D(analyticalSolutionResult, explicitSolutionResult)
	meanAbsErr, _ := meanAbsError(analyticalSolutionResult, explicitSolutionResult)
	fmt.Printf("max abs error = %.12f\n", maxAbsErr)
	fmt.Printf("mean abs error = %.12f\n", meanAbsErr)

	// Implicit
	fmt.Println("\nImplicit")
	implicitSolutionResult, _ := implicitFiniteDifferenceMethod(xRange, tRange, h, sigma, a)
	maxAbsErr, _ = maxAbsError2D(analyticalSolutionResult, implicitSolutionResult)
	meanAbsErr, _ = meanAbsError(analyticalSolutionResult, implicitSolutionResult)
	fmt.Printf("max abs error = %.12f\n", maxAbsErr)
	fmt.Printf("mean abs error = %.12f\n", meanAbsErr)

	// Crank Nicolson
	fmt.Println("\nCrank Nicolson")
	crankNicolsonResult, _ := crankNicolsonMethod(xRange, tRange, h, sigma, a, 0.5)
	maxAbsErr, _ = maxAbsError2D(analyticalSolutionResult, crankNicolsonResult)
	meanAbsErr, _ = meanAbsError(analyticalSolutionResult, crankNicolsonResult)
	fmt.Printf("max abs error = %.12f\n", maxAbsErr)
	fmt.Printf("mean abs error = %.12f\n", meanAbsErr)

	// Plot
	plotResults(map[string]struct {
		arr   [][]float64
		color color.Color
	}{
		"analytical":     {analyticalSolutionResult, color.RGBA{255, 0, 0, 255}},
		"explicit":       {explicitSolutionResult, color.RGBA{0, 255, 0, 255}},
		"implicit":       {implicitSolutionResult, color.RGBA{0, 0, 255, 255}},
		"crank-nicolson": {crankNicolsonResult, color.RGBA{255, 0, 255, 255}},
	}, 0.5, xRange, tRange, h, sigma, a)

	plotErrorsFromTime(map[string]struct {
		arr   [][]float64
		color color.Color
	}{
		"analytical":     {analyticalSolutionResult, color.RGBA{255, 0, 0, 255}},
		"explicit":       {explicitSolutionResult, color.RGBA{0, 255, 0, 255}},
		"implicit":       {implicitSolutionResult, color.RGBA{0, 0, 255, 255}},
		"crank-nicolson": {crankNicolsonResult, color.RGBA{255, 0, 255, 255}},
	}, tRange, h, sigma, a)
}
