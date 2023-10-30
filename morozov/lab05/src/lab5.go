package main

import (
	"encoding/csv"
	"fmt"
	"math"
	"os"
)

type ParabolicSolver struct {
	L                  float64
	Psi                func(float64) float64
	Function           func(float64, float64) float64
	Phi0               func(float64) float64
	PhiL               func(float64) float64
	AnaliticalSolution func(float64, float64) float64
	Type               string
}

func (solver *ParabolicSolver) Implicit(N, K int, T float64) [][]float64 {
	h := solver.L / float64(N)
	tau := T / float64(K)
	sigma := tau / (h * h)
	u := make([][]float64, K)
	for i := 0; i < K; i++ {
		u[i] = make([]float64, N)
	}

	a := make([]float64, N)
	b := make([]float64, N)
	c := make([]float64, N)
	d := make([]float64, N)

	for i := 1; i < N-1; i++ {
		u[0][i] = solver.Psi(float64(i) * h)
	}
	u[0][N-1] = 0

	for k := 1; k < K; k++ {
		for j := 1; j < N-1; j++ {
			a[j] = sigma
			b[j] = -(1 + 2*sigma)
			c[j] = sigma
			d[j] = -u[k-1][j] - tau*solver.Function(float64(j)*h, float64(k)*tau)
		}

		a[0] = 0
		b[0] = -(1 + 2*sigma)
		c[0] = sigma
		a[N-1] = sigma
		b[N-1] = -(1 + 2*sigma)
		c[N-1] = 0

		if solver.Type == "1-1" {
			d[0] = -(u[k-1][0] + sigma*solver.Phi0(float64(k)*tau))
			d[N-1] = -(u[k-1][N-1] + sigma*solver.PhiL(float64(k)*tau))
		} else if solver.Type == "2-2" {
			d[0] = -(u[k-1][0] + sigma*solver.Phi0(float64(k)*tau) - tau*solver.Function(0, float64(k)*tau))
			d[N-1] = -(u[k-1][N-1] + sigma*solver.PhiL(float64(k)*tau) - tau*solver.Function(float64(N-1)*h, float64(k)*tau))
		} else if solver.Type == "2-3" {
			d[0] = -((1-sigma)*u[k-1][1] + sigma/2*u[k-1][0] - tau*solver.Function(0, float64(k)*tau) - sigma*solver.Phi0(float64(k)*tau))
			d[N-1] = solver.PhiL(float64(k) * tau)
			d[N-1] += (tau * solver.Function(float64(N-1)*h, float64(k)*tau) * h) / (2*tau*u[k-1][N-1] + 1)
		}

		for i := 1; i < N; i++ {
			m := a[i] / b[i-1]
			b[i] -= m * c[i-1]
			d[i] -= m * d[i-1]
		}

		u[k][N-1] = d[N-1] / b[N-1]
		for i := N - 2; i >= 0; i-- {
			u[k][i] = (d[i] - c[i]*u[k][i+1]) / b[i]
		}
	}

	return u
}

func (solver *ParabolicSolver) Explicit(N, K int, T float64) [][]float64 {
	h := solver.L / float64(N)
	tau := T / float64(K)
	sigma := tau / (h * h)
	u := make([][]float64, K)
	for i := 0; i < K; i++ {
		u[i] = make([]float64, N)
	}

	for j := 1; j < N-1; j++ {
		u[0][j] = solver.Psi(float64(j) * h)
	}

	for k := 1; k < K; k++ {
		u[k][0] = solver.Phi0(float64(k) * tau)
		for j := 1; j < N-1; j++ {
			u[k][j] = sigma*u[k-1][j+1] + (1-2*sigma)*u[k-1][j] + sigma*u[k-1][j-1] + tau*solver.Function(float64(j)*h, float64(k)*tau)
		}

		if solver.Type == "1-1" {
			u[k][N-1] = u[k][N-2] + solver.PhiL(float64(k)*tau)*h
		} else if solver.Type == "2-2" {
			u[k][N-1] = solver.PhiL(float64(k) * tau)
		} else if solver.Type == "2-3" {
			u[k][N-1] = (solver.PhiL(float64(k)*tau) + u[k][N-2]/h + 2*tau*u[k-1][N-1]/h) / (1/h + 2*tau/h)
		}
	}

	return u
}

func (solver *ParabolicSolver) CrankNicholson(N, K int, T float64) [][]float64 {
	h := solver.L / float64(N)
	tau := T / float64(K)
	sigma := tau / (h * h)
	u := make([][]float64, K)
	for i := 0; i < K; i++ {
		u[i] = make([]float64, N)
	}

	for j := 1; j < N-1; j++ {
		u[0][j] = solver.Psi(float64(j) * h)
	}

	for k := 1; k < K; k++ {
		a := make([]float64, N)
		b := make([]float64, N)
		c := make([]float64, N)
		d := make([]float64, N)
		uImplicit := make([]float64, N)

		for j := 1; j < N-1; j++ {
			a[j] = sigma
			b[j] = -(1 + 2*sigma)
			c[j] = sigma
			d[j] = -u[k-1][j] - tau*solver.Function(float64(j)*h, float64(k)*tau)
		}

		a[0] = 0
		b[0] = -(1 + 2*sigma)
		c[0] = sigma
		a[N-1] = sigma
		b[N-1] = -(1 + 2*sigma)
		c[N-1] = 0

		if solver.Type == "1-1" {
			d[0] = -(u[k-1][0] + sigma*solver.Phi0(float64(k)*tau))
			d[N-1] = -(u[k-1][N-1] + sigma*solver.PhiL(float64(k)*tau))
		} else if solver.Type == "2-2" {
			d[0] = -(u[k-1][0] + sigma*solver.Phi0(float64(k)*tau) - tau*solver.Function(0, float64(k)*tau))
			d[N-1] = -(u[k-1][N-1] + sigma*solver.PhiL(float64(k)*tau) - tau*solver.Function(float64(N-1)*h, float64(k)*tau))
		} else if solver.Type == "2-3" {
			d[0] = -((1-sigma)*u[k-1][1] + sigma/2*u[k-1][0] - tau*solver.Function(0, float64(k)*tau) - sigma*solver.Phi0(float64(k)*tau))
			d[N-1] = solver.PhiL(float64(k) * tau)
			d[N-1] += (tau * solver.Function(float64(N-1)*h, float64(k)*tau) * h) / (2*tau*u[k-1][N-1] + 1)
		}

		p := make([]float64, N)
		q := make([]float64, N)
		p[0] = -c[0] / b[0]
		q[0] = d[0] / b[0]

		for i := 1; i < N; i++ {
			m := a[i] / b[i-1]
			b[i] -= m * c[i-1]
			d[i] -= m * d[i-1]
		}

		uImplicit[N-1] = q[N-1]
		for i := N - 2; i >= 0; i-- {
			uImplicit[i] = p[i]*uImplicit[i+1] + q[i]
		}

		uExplicit := make([]float64, N)
		for j := 1; j < N-1; j++ {
			uExplicit[j] = sigma*u[k-1][j+1] + (1-2*sigma)*u[k-1][j] + sigma*u[k-1][j-1] + tau*solver.Function(float64(j)*h, float64(k)*tau)
		}

		if solver.Type == "1-1" {
			uImplicit[N-1] = u[k][N-2] + solver.PhiL(float64(k)*tau)*h
		} else if solver.Type == "2-2" {
			uImplicit[N-1] = solver.PhiL(float64(k) * tau)
		} else if solver.Type == "2-3" {
			uImplicit[N-1] = (solver.PhiL(float64(k)*tau) + u[k][N-2]/h + 2*tau*u[k-1][N-1]/h) / (1/h + 2*tau/h)
		}

		for j := 0; j < N; j++ {
			u[k][j] = 0.5*uImplicit[j] + 0.5*uExplicit[j]
		}
	}

	return u
}

func AnalyticalSolutionMatrix(N, K int, T float64, solver *ParabolicSolver) [][]float64 {
	h := solver.L / float64(N)
	tau := T / float64(K)
	u := make([][]float64, K)
	for k := 0; k < K; k++ {
		u[k] = make([]float64, N)
		for j := 0; j < N; j++ {
			u[k][j] = solver.AnaliticalSolution(float64(j)*h, float64(k)*tau)
		}
	}
	return u
}

func writeCSV(filename string, data [][]float64) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()
	writer := csv.NewWriter(file)
	defer writer.Flush()
	for _, row := range data {
		stringRow := make([]string, len(row))
		for i, val := range row {
			stringRow[i] = fmt.Sprintf("%f", val)
		}
		if err := writer.Write(stringRow); err != nil {
			return err
		}
	}
	return nil
}

func main() {
	solver := &ParabolicSolver{
		L:                  math.Pi,
		Psi:                func(x float64) float64 { return math.Sin(x) },
		Function:           func(x, t float64) float64 { return 0.5 * math.Exp(-0.5*t) * math.Cos(x) },
		Phi0:               func(t float64) float64 { return -math.Exp(-0.5 * t) },
		PhiL:               func(t float64) float64 { return -math.Exp(-0.5 * t) },
		AnaliticalSolution: func(x, t float64) float64 { return math.Exp(-0.5*t) * math.Sin(x) },
		Type:               "1-1",
	}
	T, K, N := 10.0, 500, 5
	algorithms := []string{"Implicit", "Explicit", "CrankNicholson", "Analytic"}
	answers := make(map[string][][]float64)
	for _, algorithm := range algorithms {
		if algorithm == "Implicit" {
			answers[algorithm] = solver.Implicit(N, K, T)
		} else if algorithm == "Explicit" {
			answers[algorithm] = solver.Explicit(N, K, T)
		} else if algorithm == "CrankNicholson" {
			answers[algorithm] = solver.CrankNicholson(N, K, T)
		} else if algorithm == "Analytic" {
			answers[algorithm] = AnalyticalSolutionMatrix(N, K, T, solver)
		}
	}
	for _, algorithm := range algorithms {
		data := answers[algorithm]
		filename := algorithm + ".csv"
		err := writeCSV(filename, data)
		if err != nil {
			panic(err)
		}
	}
}
