package main

import "math"

type Interval struct {
	Start, End float64
}

func (interval *Interval) Range(step float64) []float64 {
	size := int(math.Ceil((math.Max(interval.Start, interval.End) - math.Min(interval.Start, interval.End)) / step))
	rng := make([]float64, size)
	for i := range rng {
		rng[i] = math.Min(interval.Start, interval.End) + float64(i)*step
	}
	return rng
}

func transpose3D(U [][][]float64) [][][]float64 {
	if len(U) == 0 || len(U[0]) == 0 || len(U[0][0]) == 0 {
		return U // returns the original slice if it's empty or has zero length in any dimension
	}

	K := len(U)
	Nx := len(U[0])
	Ny := len(U[0][0])
	UT := make([][][]float64, Ny) // Transposed slice

	for i := 0; i < Ny; i++ {
		UT[i] = make([][]float64, Nx)
		for j := 0; j < Nx; j++ {
			UT[i][j] = make([]float64, K)
			for k := 0; k < K; k++ {
				UT[i][j][k] = U[k][j][i]
			}
		}
	}

	return UT
}
