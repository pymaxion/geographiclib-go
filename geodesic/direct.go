package geodesic

import (
	"math"

	"geographiclib-go/geodesic/capabilities"
)

type directSolver struct {
	g *geodesicImpl
}

func (s *directSolver) direct(lat1, lon1, azi1, s12 float64, mask capabilities.BitMask) *dataImpl {
	return s.genDirect(lat1, lon1, azi1, false, s12, mask)
}

func (s *directSolver) arcDirect(lat1, lon1, azi1, a12 float64, mask capabilities.BitMask) *dataImpl {
	return s.genDirect(lat1, lon1, azi1, true, a12, mask)
}

//goland:noinspection GoSnakeCaseUsage
func (s *directSolver) genDirect(lat1, lon1, azi1 float64, arcMode bool, s12_a12 float64, mask capabilities.BitMask) *dataImpl {
	// Automatically supply DistanceIn if necessary
	if !arcMode {
		mask |= capabilities.DistanceIn
	}
	line := newLineImpl(s.g, lat1, lon1, azi1, math.NaN(), math.NaN(), mask)
	return line.solvePosition(arcMode, s12_a12, mask)
}
