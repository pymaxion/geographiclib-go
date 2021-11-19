package geodesic

import (
	"geographiclib-go/geodesic/capabilities"
)

type PolygonArea struct {
	g            *Geodesic
	polyline     bool
	area0        float64
	caps         capabilities.Mask
	perimeterSum *accumulator
	areaSum      *accumulator
}

func newPolygonArea(g *Geodesic, polyline bool) *PolygonArea {
	area0 := g.EllipsoidArea()
	caps := capabilities.Latitude | capabilities.Longitude | capabilities.Distance
	if !polyline {
		caps |= capabilities.Area | capabilities.LongUnroll
	}
	perimeterSum := newAccumulator(0)
	var areaSum *accumulator
	if !polyline {
		areaSum = newAccumulator(0)
	}

	p := &PolygonArea{
		g:            g,
		polyline:     polyline,
		area0:        area0,
		caps:         caps,
		perimeterSum: perimeterSum,
		areaSum:      areaSum,
	}

	p.Clear()
	return p
}

func (p *PolygonArea) Clear() {
	//p.num = 0
	//p.crossings = 0
	//p.perimeterSum.set(0)
	//if !p.polyline {
	//	p.areaSum.set(0)
	//}
	//p.lat0, p.lon0, p.lat1, p.lon1 = math.NaN(), math.NaN(), math.NaN(), math.NaN()
}
