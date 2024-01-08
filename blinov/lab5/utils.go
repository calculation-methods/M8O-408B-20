package main

import (
	"fmt"
	"math"
)

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

func maxAbsError2D(A, B [][]float64) (float64, error) {
	if len(A) == 0 || len(B) == 0 || len(A) != len(B) || len(A[0]) != len(B[0]) {
		return 0, fmt.Errorf("matrices must be non-empty and have the same dimensions")
	}

	var maxError float64 = 0.0
	for i := range A {
		if len(A[i]) != len(B[i]) {
			return 0, fmt.Errorf("matrices must have the same dimensions")
		}
		for j := range A[i] {
			vError := math.Abs(A[i][j] - B[i][j])
			if vError > maxError {
				maxError = vError
			}
		}
	}

	return maxError, nil
}

func maxAbsError1D(A, B []float64) (float64, error) {
	if len(A) == 0 || len(B) == 0 || len(A) != len(B) {
		return 0, fmt.Errorf("matrices must be non-empty and have the same dimensions")
	}

	var maxError float64 = 0.0
	for i := range A {
		vError := math.Abs(A[i] - B[i])
		if vError > maxError {
			maxError = vError
		}
	}

	return maxError, nil
}

func meanAbsError(A, B [][]float64) (float64, error) {
	if len(A) == 0 || len(B) == 0 || len(A) != len(B) || len(A[0]) != len(B[0]) {
		return 0, fmt.Errorf("matrices must be non-empty and have the same dimensions")
	}

	totalError := 0.0
	count := 0

	for i := range A {
		if len(A[i]) != len(B[i]) {
			return 0, fmt.Errorf("matrices must have the same dimensions")
		}
		for j := range A[i] {
			vError := math.Abs(A[i][j] - B[i][j])
			totalError += vError
			count++
		}
	}

	if count == 0 {
		return 0, fmt.Errorf("no elements to compare")
	}

	meanError := totalError / float64(count)
	return meanError, nil
}
