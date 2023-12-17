package main

import (
	"fmt"
	"image/color"
	"math"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
)

func plotResults(solutions map[string]struct {
	arr   [][]float64
	color color.Color
}, time float64, xRange, tRange Interval, h, sigma, a float64) error {
	tau := sigma * math.Pow(h, 2) / a // len of cell by t
	x := xRange.Range(h)
	times := tRange.Range(tau)
	curTID := findClosestIndex(times, time)

	p := plot.New()
	p.Title.Text = "Solutions"
	p.X.Label.Text = "X"
	p.Y.Label.Text = "T"
	p.Add(plotter.NewGrid())

	for methodName, solution := range solutions {
		pts := make(plotter.XYs, len(x))
		for i, xVal := range x {
			pts[i].X = xVal
			pts[i].Y = solution.arr[curTID][i]
		}
		line, err := plotter.NewLine(pts)
		if err != nil {
			return fmt.Errorf("could not create line for %s: %v", methodName, err)
		}
		line.Color = solution.color
		p.Legend.Add(methodName, line)
		line.Width = 3
		p.Add(line)
	}

	p.Legend.Top = true
	p.Legend.TextStyle.Font.Size = 18
	if err := p.Save(12*vg.Inch, 6*vg.Inch, "solutions.png"); err != nil {
		return fmt.Errorf("could not save plot: %v", err)
	}
	return nil
}

func findClosestIndex(times []float64, time float64) int {
	minDiff := math.MaxFloat64
	minIndex := 0
	for i, t := range times {
		diff := math.Abs(t - time)
		if diff < minDiff {
			minDiff = diff
			minIndex = i
		}
	}
	return minIndex
}

func plotErrorsFromTime(
	solutions map[string]struct {
		arr   [][]float64
		color color.Color
	},
	tRange Interval,
	h, sigma, a float64,
) error {
	tau := sigma * math.Pow(h, 2) / a
	t := tRange.Range(tau)

	p := plot.New()
	p.Title.Text = "Max Absolute Error over Time"
	p.X.Label.Text = "Time"
	p.Y.Label.Text = "Max Absolute Error"
	p.Add(plotter.NewGrid())

	for methodName, solution := range solutions {
		if methodName == "analytical" {
			continue
		}

		pts := make(plotter.XYs, len(t))
		for i, time := range t {
			pts[i].X = time
			pts[i].Y, _ = maxAbsError1D(solution.arr[i], solutions["analytical"].arr[i])
		}

		line, err := plotter.NewLine(pts)
		if err != nil {
			return fmt.Errorf("could not create line for %s: %v", methodName, err)
		}
		line.Color = solution.color
		line.Width = 3
		p.Legend.Add(methodName, line)
		p.Add(line)
	}

	p.Legend.Top = true

	if err := p.Save(12*vg.Inch, 6*vg.Inch, "error_plot.png"); err != nil {
		return fmt.Errorf("could not save plot: %v", err)
	}
	return nil
}
