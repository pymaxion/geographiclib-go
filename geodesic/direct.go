package geodesic

import (
	"math"

	"github.com/pymaxion/geographiclib-go/geodesic/capabilities"
)

type directSolver struct {
	g *Geodesic
}

func (s *directSolver) direct(lat1, lon1, azi1, s12 float64, caps capabilities.Mask) Data {
	caps |= capabilities.DistanceIn // automatically supply DistanceIn if necessary
	line := newLine(s.g, lat1, lon1, azi1, math.NaN(), math.NaN(), caps)
	return line.solvePosition(false, s12, caps)
}

func (s *directSolver) arcDirect(lat1, lon1, azi1, a12 float64, caps capabilities.Mask) Data {
	line := newLine(s.g, lat1, lon1, azi1, math.NaN(), math.NaN(), caps)
	return line.solvePosition(true, a12, caps)
}
