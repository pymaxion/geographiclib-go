package geodesic

import (
	"math"

	"geographiclib-go/geodesic/capabilities"
)

type PolygonArea struct {
	g            *Geodesic
	polyline     bool
	area0        float64
	caps         capabilities.Mask
	perimeterSum *accumulator
	areaSum      *accumulator
	num          int
	crossings    int
	lat0         float64
	lon0         float64
	lat1         float64
	lon1         float64
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
	p.num = 0
	p.crossings = 0
	p.perimeterSum.set(0)
	if !p.polyline {
		p.areaSum.set(0)
	}
	p.lat0, p.lon0, p.lat1, p.lon1 = math.NaN(), math.NaN(), math.NaN(), math.NaN()
}

func (p *PolygonArea) AddPoint(lat, lon float64) {
	lon = angNormalize(lon)
	if p.num == 0 {
		p.lat0, p.lat1 = lat, lat
		p.lon0, p.lon1 = lon, lon
	} else {
		inv := p.g.InverseWithCapabilities(p.lat1, p.lon1, lat, lon, p.caps)
		p.perimeterSum.add(inv.S12)
		if !p.polyline {
			p.areaSum.add(inv.S12Area)
			p.crossings += transit(p.lon1, lon)
		}
		p.lat1, p.lon1 = lat, lon
	}
	p.num++
}

func (p *PolygonArea) AddEdge(azi, s float64) {
	if p.num > 0 { // Do nothing if p.num is zero
		dir := p.g.DirectWithCapabilities(p.lat1, p.lon1, azi, s, p.caps)
		p.perimeterSum.add(dir.S12)
		if !p.polyline {
			p.areaSum.add(dir.S12Area)
			p.crossings += transitDirect(p.lon1, dir.Lon2)
		}
		p.lat1, p.lon1 = dir.Lat2, dir.Lon2
		p.num++
	}
}

type PolygonResult struct {
	num       int
	perimeter float64
	area      float64
}

func (p *PolygonArea) Compute(reverse, sign bool) PolygonResult {
	if p.num < 2 {
		return PolygonResult{p.num, 0, ternary(p.polyline, math.NaN(), 0)}
	}
	if p.polyline {
		return PolygonResult{p.num, p.perimeterSum.sum(), math.NaN()}
	}

	inv := p.g.InverseWithCapabilities(p.lat1, p.lon1, p.lat0, p.lon0, p.caps)
	tempSum := *p.areaSum
	tempSum.add(inv.S12Area)

	area := areaReduceA(&tempSum, p.area0, p.crossings+transit(p.lon1, p.lon0), reverse, sign)
	return PolygonResult{p.num, p.perimeterSum.sumWith(inv.S12), area}
}

func (p *PolygonArea) TestPoint(lat, lon float64, reverse, sign bool) PolygonResult {
	if p.num == 0 {
		return PolygonResult{1, 0, ternary(p.polyline, math.NaN(), 0)}
	}

	perimeter := p.perimeterSum.sum()
	tempSum := 0.
	if !p.polyline {
		tempSum = p.areaSum.sum()
	}
	crossings := p.crossings
	num := p.num + 1
	endIdx := int(ternary(p.polyline, 1, 2))
	for i := 0; i < endIdx; i++ {
		inv := p.g.InverseWithCapabilities(
			ternary(i == 0, p.lat1, lat),
			ternary(i == 0, p.lon1, lon),
			ternary(i != 0, p.lat0, lat),
			ternary(i != 0, p.lon0, lon),
			p.caps)
		perimeter += inv.S12
		if !p.polyline {
			tempSum += inv.S12Area
			crossings += transit(ternary(i == 0, p.lon1, lon), ternary(i != 0, p.lon0, lon))
		}
	}

	if p.polyline {
		return PolygonResult{num, perimeter, math.NaN()}
	}

	area := areaReduceB(tempSum, p.area0, crossings, reverse, sign)
	return PolygonResult{num, perimeter, area}
}

func (p *PolygonArea) TestEdge(azi, s float64, reverse, sign bool) PolygonResult {
	if p.num == 0 { // we don't have a starting point!
		return PolygonResult{0, math.NaN(), math.NaN()}
	}

	num := p.num + 1
	perimeter := p.perimeterSum.sum() + s
	if p.polyline {
		return PolygonResult{num, perimeter, math.NaN()}
	}

	tempSum := p.areaSum.sum()
	crossings := p.crossings
	dir := p.g.DirectWithCapabilities(p.lat1, p.lon1, azi, s, p.caps)
	tempSum += dir.S12Area
	crossings += transitDirect(p.lon1, dir.Lon2)
	crossings += transit(dir.Lon2, p.lon0)
	inv := p.g.InverseWithCapabilities(dir.Lat2, dir.Lon2, p.lat0, p.lon0, p.caps)
	perimeter += inv.S12
	tempSum += inv.S12Area

	area := areaReduceB(tempSum, p.area0, crossings, reverse, sign)
	return PolygonResult{num, perimeter, area}
}

// transit returns 1 or -1 if crossing prime meridian in east or west direction, else zero.
func transit(lon1, lon2 float64) int {
	// Compute lon12 the same way as Geodesic.Inverse.
	lon1 = angNormalize(lon1)
	lon2 = angNormalize(lon2)
	lon12, _ := angDiff(lon1, lon2)
	cross := 0
	if lon1 <= 0 && lon2 > 0 && lon12 > 0 {
		cross = 1
	} else if lon2 <= 0 && lon1 > 0 && lon12 < 0 {
		cross = -1
	}
	return cross
}

// transitDirect is an alternate version of transit to deal with longitudes in the direct problem.
func transitDirect(lon1, lon2 float64) int {
	// We want to compute exactly:
	//   int(ceil(lon2 / 360)) - int(ceil(lon1 / 360))
	// Since we only need the parity of the result we can use std::remquo but this is buggy with g++
	// 4.8.3 and requires C++11. So instead we do
	lon1 = math.Mod(lon1, 720.0)
	lon2 = math.Mod(lon2, 720.0)
	u, v := 0, 0
	if (lon2 <= 0 && lon2 > -360) || lon2 > 360 {
		u = 1
	}
	if (lon1 <= 0 && lon1 > -360) || lon1 > 360 {
		v = 1
	}
	return u - v
}

func areaReduceA(area *accumulator, area0 float64, crossings int, reverse, sign bool) float64 {
	area.remainder(area0)
	if (crossings & 1) != 0 {
		area.add(ternary(area.sum() < 0, 1, -1) * area0 / 2)
	}
	// area is with the clockwise sense. If !reverse convert to counter-clockwise convention.
	if !reverse {
		area.negate()
	}
	// If sign put area in (-area0/2, area0/2], else put area in [0, area0)
	if sign {
		if area.sum() > area0/2 {
			area.add(-area0)
		} else if area.sum() <= -area0/2 {
			area.add(+area0)
		}
	} else {
		if area.sum() >= area0 {
			area.add(-area0)
		} else if area.sum() < 0 {
			area.add(+area0)
		}
	}
	return 0 + area.sum()
}

func areaReduceB(area, area0 float64, crossings int, reverse, sign bool) float64 {
	area = remainder(area, area0)
	if (crossings & 1) != 0 {
		area += ternary(area < 0, 1, -1) * area0 / 2
	}
	// area is with the clockwise sense. If !reverse convert to counter-clockwise convention.
	if !reverse {
		area *= -1
	}
	// If sign put area in (-area0/2, area0/2], else put area in [0, area0)
	if sign {
		if area > area0/2 {
			area -= area0
		} else if area <= -area0/2 {
			area += area0
		}
	} else {
		if area >= area0 {
			area -= area0
		} else if area < 0 {
			area += area0
		}
	}
	return 0 + area
}
