package geodesic

import (
	"math"

	"github.com/pymaxion/geographiclib-go/geodesic/capabilities"
)

type Gnomonic struct {
	Earth *Geodesic
}

func (g *Gnomonic) Forward(lat0, lon0, lat, lon float64) GnomonicData {
	fwd := newGnomonicData()
	fwd.Lat0, fwd.Lon0, fwd.Lat, fwd.Lon = lat0, lon0, lat, lon

	caps := capabilities.Azimuth | capabilities.GeodesicScale | capabilities.ReducedLength
	inv := g.Earth.InverseWithCapabilities(lat0, lon0, lat, lon, caps)
	fwd.Azi, fwd.Rk = inv.Azi2, inv.M12

	if inv.M12 > 0 {
		rho := inv.M12Reduced / inv.M12
		x, y := sincosd(inv.Azi1)
		fwd.X, fwd.Y = rho*x, rho*y
	}
	return fwd
}

func (g *Gnomonic) Reverse(lat0, lon0, x, y float64) GnomonicData {
	rev := newGnomonicData()
	rev.Lat0, rev.Lon0, rev.X, rev.Y = lat0, lon0, x, y

	azi0 := atan2d(x, y)
	rho := math.Hypot(x, y)
	a := g.Earth.EquatorialRadius()
	s := a * math.Atan(rho/a)

	little := rho <= a
	if !little {
		rho = 1 / rho
	}

	caps := capabilities.Latitude | capabilities.Longitude | capabilities.Azimuth | capabilities.DistanceIn |
		capabilities.ReducedLength | capabilities.GeodesicScale
	line := g.Earth.LineWithCapabilities(lat0, lon0, azi0, caps)

	var pos Data
	numIt := 10
	trip := 0
	tripEpsilon := 0.01 * math.Sqrt(epsilon)

	for count := numIt; count > 0; count-- {
		pos = line.PositionWithCapabilities(s, caps)
		if trip > 0 {
			break
		}

		var ds float64
		if little {
			ds = ((pos.M12Reduced / pos.M12) - rho) * pos.M12 * pos.M12
		} else {
			ds = (rho - (pos.M12 / pos.M12Reduced)) * pos.M12Reduced * pos.M12Reduced
		}
		s -= ds

		if math.Abs(ds) < tripEpsilon*a {
			trip++
		}
	}

	if trip == 0 {
		return rev
	}

	rev.Lat, rev.Lon, rev.Azi, rev.Rk = pos.Lat2, pos.Lon2, pos.Azi1, pos.M12
	return rev
}

type GnomonicData struct {
	Lat0 float64
	Lon0 float64
	Lat  float64
	Lon  float64
	X    float64
	Y    float64
	Azi  float64
	Rk   float64
}

func newGnomonicData() GnomonicData {
	return GnomonicData{
		Lat0: math.NaN(),
		Lon0: math.NaN(),
		Lat:  math.NaN(),
		Lon:  math.NaN(),
		X:    math.NaN(),
		Y:    math.NaN(),
		Azi:  math.NaN(),
		Rk:   math.NaN(),
	}
}
