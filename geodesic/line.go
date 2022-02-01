package geodesic

import (
	"math"

	"geographiclib-go/geodesic/capabilities"
)

/*
 Line represents a geodesic line and facilitates the determination of a series of points on a single
 geodesic. Geodesic.Line method should be used to create an instance of Line.

 Position returns the location of point 2 a distance s12 along the geodesic. Alternatively,
 ArcPosition gives the position of point 2 an arc length a12 along the geodesic. The additional
 functions PositionWithCapabilities and ArcPositionWithCapabilities include an optional final
 argument of type capabilities.Mask to allow you to specify which results should be computed and
 returned.

 You can register the position of a reference point 3 a distance (arc length), s13 (a13) along the
 geodesic with the SetDistance (SetArc) functions. Points a fractional distance along the line can
 be found by providing, for example, 0.5 * Distance as an argument to Position. The
 Geodesic.InverseLine or Geodesic.DirectLine functions return Line instances with point 3 set to the
 point 2 of the corresponding geodesic problem. Line instances created Geodesic.Line have s13 and
 a13 set to math.NaN.

 The calculations are accurate to better than 15 nm (15 nanometers). See Sec. 9 of arXiv:1102.1215v1
 (https://arxiv.org/abs/1102.1215v1) for details. The algorithms used by this class are based on
 series expansions using the flattening f as a small parameter. These are only accurate for |f| <
 0.02; however reasonably accurate results will be obtained for |f| < 0.2.

 The algorithms are described in
  C. F. F. Karney, Algorithms for geodesics, J. Geodesy 87, 43-55 (2013)
  Link: https://doi.org/10.1007/s00190-012-0578-z
  Addenda: https://geographiclib.sourceforge.io/geod-addenda.html
*/
type Line struct {
	g     *Geodesic
	lat1  float64
	lon1  float64
	azi1  float64
	salp1 float64
	calp1 float64
	dn1   float64
	salp0 float64
	calp0 float64
	ssig1 float64
	csig1 float64
	somg1 float64
	comg1 float64
	k2    float64
	stau1 float64
	ctau1 float64
	a1m1  float64
	a2m1  float64
	a3c   float64
	b11   float64
	b21   float64
	b31   float64
	a4    float64
	b41   float64
	a13   float64
	s13   float64
	c1a   []float64
	c1pa  []float64
	c2a   []float64
	c3a   []float64
	c4a   []float64
	mask  capabilities.Mask
}

func newLine(g *Geodesic, lat1, lon1, azi1, salp1, calp1 float64, caps capabilities.Mask) *Line {
	// Always allow latitude and azimuth and unrolling the longitude
	caps |= capabilities.Latitude | capabilities.Azimuth | capabilities.LongUnroll
	lat1 = latFix(lat1)
	if math.IsNaN(salp1) || math.IsNaN(calp1) {
		azi1 = angNormalize(azi1)
		salp1, calp1 = sincosd(angRound(azi1))
	}

	sbet1, cbet1 := sincosd(angRound(lat1))
	sbet1 *= g.f1
	// Ensure cbet1 = +epsilon at poles
	sbet1, cbet1 = norm(sbet1, cbet1)
	cbet1 = math.Max(tiny, cbet1)
	dn1 := math.Sqrt(1 + g.ep2*sq(sbet1))

	// Evaluate alp0 from sin(alp1) * cos(bet1) = sin(alp0),
	salp0 := salp1 * cbet1 // alp0 in [0, pi/2 - |bet1|]
	// Alt: calp0 = Math.hypot(sbet1, calp1 * cbet1).  The following
	// is slightly better (consider the case salp1 = 0).
	calp0 := math.Hypot(calp1, salp1*sbet1)
	// Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
	// sig = 0 is nearest northward crossing of equator.
	// With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
	// With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
	// With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
	// Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
	// With alp0 in (0, pi/2], quadrants for sig and omg coincide.
	// No atan2(0,0) ambiguity at poles since cbet1 = +epsilon.
	// With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi.
	ssig1 := sbet1
	somg1 := salp0 * sbet1
	csig1 := 1.
	if sbet1 != 0 || calp1 != 0 {
		csig1 = cbet1 * calp1
	}
	comg1 := csig1
	ssig1, csig1 = norm(ssig1, csig1)
	// somg1, comg1 = norm(somg1, comg1) -- don't need to normalize these!

	k2 := sq(calp0) * g.ep2
	eps := k2 / (2*(1+math.Sqrt(1+k2)) + k2)

	var stau1, ctau1, a1m1, a2m1, a3c, b11, b21, b31, a4, b41 = math.NaN(), math.NaN(), math.NaN(),
		math.NaN(), math.NaN(), math.NaN(), math.NaN(), math.NaN(), math.NaN(), math.NaN()
	var c1a, c1pa, c2a, c3a, c4a []float64

	if (caps & capabilities.C1) != 0 {
		a1m1 = a1m1f(eps)
		c1a = make([]float64, nC1+1)
		c1f(eps, c1a)
		b11 = sinCosSeries(true, ssig1, csig1, c1a)
		s, c := math.Sincos(b11)
		// tau1 = sig1 + B11
		stau1 = ssig1*c + csig1*s
		ctau1 = csig1*c - ssig1*s
		// Not necessary because C1pa reverts C1a
		// b11 = -sinCosSeries(true, stau1, ctau1, c1pa, nC1p)
	}

	if (caps & capabilities.C1p) != 0 {
		c1pa = make([]float64, nC1p+1)
		c1pf(eps, c1pa)
	}

	if (caps & capabilities.C2) != 0 {
		c2a = make([]float64, nC2+1)
		a2m1 = a2m1f(eps)
		c2f(eps, c2a)
		b21 = sinCosSeries(true, ssig1, csig1, c2a)
	}

	if (caps & capabilities.C3) != 0 {
		c3a = make([]float64, nC3)
		g.c3f(eps, c3a)
		a3c = -g.f * salp0 * g.a3f(eps)
		b31 = sinCosSeries(true, ssig1, csig1, c3a)
	}

	if (caps & capabilities.C4) != 0 {
		c4a = make([]float64, nC4)
		g.c4f(eps, c4a)
		// Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0)
		a4 = sq(g.a) * calp0 * salp0 * g.e2
		b41 = sinCosSeries(false, ssig1, csig1, c4a)
	}

	return &Line{
		g:     g,
		lat1:  lat1,
		lon1:  lon1,
		azi1:  azi1,
		salp1: salp1,
		calp1: calp1,
		dn1:   dn1,
		salp0: salp0,
		calp0: calp0,
		ssig1: ssig1,
		csig1: csig1,
		somg1: somg1,
		comg1: comg1,
		k2:    k2,
		stau1: stau1,
		ctau1: ctau1,
		a1m1:  a1m1,
		a2m1:  a2m1,
		a3c:   a3c,
		b11:   b11,
		b21:   b21,
		b31:   b31,
		a4:    a4,
		b41:   b41,
		a13:   math.NaN(),
		s13:   math.NaN(),
		c1a:   c1a,
		c1pa:  c1pa,
		c2a:   c2a,
		c3a:   c3a,
		c4a:   c4a,
		mask:  caps,
	}
}

/*
 Position computes the position of point 2 which is a distance s12 (meters) from point 1. The values
 of lon2 and azi2 returned are in the range [-180°, 180°].

  s12: distance from point 1 to point 2 (meters); can be negative.

 This function is equivalent to calling PositionWithCapabilities with capabilities.Standard.
*/
func (l *Line) Position(s12 float64) Data {
	return l.PositionWithCapabilities(s12, capabilities.Standard)
}

/*
 PositionWithCapabilities computes the position of point 2 which is a distance s12 (meters) from
 point 1. It also allows you to specify which results should be computed and returned via the
 capabilities.Mask argument. Note that the Line instance must have been created with caps |=
 capabilities.DistanceIn; otherwise, no parameters are set.

 See Position for more details.
*/
func (l *Line) PositionWithCapabilities(s12 float64, caps capabilities.Mask) Data {
	return l.solvePosition(false, s12, caps)
}

/*
 ArcPosition computes the position of point 2 which is an arc length a12 (degrees) from point 1. The
 values of lon2 and azi2 returned are in the range [-180°, 180°].

  a12: arc length from point 1 to point 2 (degrees); can be negative.

 This function is equivalent to calling ArcPositionWithCapabilities with capabilities.Standard.
*/
func (l *Line) ArcPosition(a12 float64) Data {
	return l.ArcPositionWithCapabilities(a12, capabilities.Standard)
}

/*
 ArcPositionWithCapabilities computes the position of point 2 which is an arc length a12 (degrees)
 from point 1. It also allows you to specify which results should be computed and returned via the
 capabilities.Mask argument. Note that the Line instance must have been created with caps |=
 capabilities.DistanceIn; otherwise, no parameters are set.

 See ArcPosition for more details.
*/
func (l *Line) ArcPositionWithCapabilities(a12 float64, caps capabilities.Mask) Data {
	return l.solvePosition(true, a12, caps)
}

//goland:noinspection GoSnakeCaseUsage
func (l *Line) solvePosition(arcMode bool, s12_a12 float64, caps capabilities.Mask) Data {
	caps &= capabilities.OutMask

	var lon1 float64
	if (caps & capabilities.LongUnroll) != 0 {
		lon1 = l.lon1
	} else {
		lon1 = angNormalize(l.lon1)
	}

	pr := l.genPosition(arcMode, s12_a12, caps)
	return Data{
		Lat1:       latFix(l.lat1),
		Lon1:       lon1,
		Azi1:       angNormalize(l.azi1),
		Lat2:       pr.lat2,
		Lon2:       pr.lon2,
		Azi2:       pr.azi2,
		S12:        pr.s12,
		A12:        pr.a12,
		M12Reduced: pr.m12,
		M12:        pr.M12,
		M21:        pr.M21,
		S12Area:    pr.S12,
	}
}

type positionResult struct {
	a12  float64
	lat2 float64
	lon2 float64
	azi2 float64
	s12  float64
	m12  float64
	M12  float64
	M21  float64
	S12  float64
}

//goland:noinspection GoSnakeCaseUsage
func (l *Line) genPosition(arcMode bool, s12_a12 float64, caps capabilities.Mask) positionResult {
	r := positionResult{
		a12:  math.NaN(),
		lat2: math.NaN(),
		lon2: math.NaN(),
		azi2: math.NaN(),
		s12:  math.NaN(),
		m12:  math.NaN(),
		M12:  math.NaN(),
		M21:  math.NaN(),
		S12:  math.NaN(),
	}

	caps &= l.mask & capabilities.OutMask
	if !arcMode && ((l.mask & capabilities.OutMask & capabilities.DistanceIn) == 0) {
		return r // Uninitialized or impossible distance calculation requested
	}

	var sig12, ssig12, csig12, b12, ab1 float64
	if arcMode {
		// Interpret s12_a12 as spherical arc length
		r.a12 = s12_a12
		sig12 = deg2rad(s12_a12)
		ssig12, csig12 = sincosd(s12_a12)
	} else {
		// Interpret s12_a12 as distance
		r.s12 = s12_a12
		tau12 := s12_a12 / (l.g.b * (1 + l.a1m1))
		s, c := math.Sincos(tau12)
		// tau2 = tau1 + tau12
		b12 = -sinCosSeries(true, l.stau1*c+l.ctau1*s, l.ctau1*c-l.stau1*s, l.c1pa)
		sig12 = tau12 - (b12 - l.b11)
		ssig12, csig12 = math.Sincos(sig12)
		if math.Abs(l.g.f) > 0.01 {
			// Reverted distance series is inaccurate for |f| > 1/100, so correct
			// sig12 with 1 Newton iteration.  The following table shows the
			// approximate maximum error for a = WGS_a() and various f relative to
			// GeodesicExact.
			//     erri = the error in the inverse solution (nm)
			//     errd = the error in the direct solution (series only) (nm)
			//     errda = the error in the direct solution
			//             (series + 1 Newton) (nm)
			//
			//       f     erri  errd errda
			//     -1/5    12e6 1.2e9  69e6
			//     -1/10  123e3  12e6 765e3
			//     -1/20   1110 108e3  7155
			//     -1/50  18.63 200.9 27.12
			//     -1/100 18.63 23.78 23.37
			//     -1/150 18.63 21.05 20.26
			//      1/150 22.35 24.73 25.83
			//      1/100 22.35 25.03 25.31
			//      1/50  29.80 231.9 30.44
			//      1/20   5376 146e3  10e3
			//      1/10  829e3  22e6 1.5e6
			//      1/5   157e6 3.8e9 280e6
			ssig2 := l.ssig1*csig12 + l.csig1*ssig12
			csig2 := l.csig1*csig12 - l.ssig1*ssig12
			b12 = sinCosSeries(true, ssig2, csig2, l.c1a)
			serr := (1+l.a1m1)*(sig12+(b12-l.b11)) - s12_a12/l.g.b
			sig12 = sig12 - serr/math.Sqrt(1+l.k2*sq(ssig2))
			ssig12, csig12 = math.Sincos(sig12)
			// Update B12 below
		}
		r.a12 = rad2deg(sig12)
	}

	var ssig2, csig2, sbet2, cbet2, salp2, calp2 float64
	// sig2 = sig1 + sig12
	ssig2 = l.ssig1*csig12 + l.csig1*ssig12
	csig2 = l.csig1*csig12 - l.ssig1*ssig12
	dn2 := math.Sqrt(1 + l.k2*sq(ssig2))
	if (caps & (capabilities.Distance | capabilities.ReducedLength | capabilities.GeodesicScale)) != 0 {
		if arcMode || math.Abs(l.g.f) > 0.01 {
			b12 = sinCosSeries(true, ssig2, csig2, l.c1a)
		}
		ab1 = (1 + l.a1m1) * (b12 - l.b11)
	}
	// sin(bet2) = cos(alp0) * sin(sig2)
	sbet2 = l.calp0 * ssig2
	// Alt: cbet2 = Math.hypot(csig2, salp0 * ssig2);
	cbet2 = math.Hypot(l.salp0, l.calp0*csig2)
	if cbet2 == 0 {
		// I.e., salp0 = 0, csig2 = 0.  Break the degeneracy in this case
		cbet2, csig2 = tiny, tiny
	}
	// tan(alp0) = cos(sig2)*tan(alp2)
	salp2, calp2 = l.salp0, l.calp0*csig2 // No need to normalize

	if ((caps & capabilities.Distance) != 0) && arcMode {
		r.s12 = l.g.b * ((1+l.a1m1)*sig12 + ab1)
	}

	if (caps & capabilities.Longitude) != 0 {
		// tan(omg2) = sin(alp0) * tan(sig2)
		somg2, comg2 := l.salp0*ssig2, csig2       // No need to normalize
		E := ternary(math.Signbit(l.salp0), -1, 1) // east or west going?
		// omg12 = omg2 - omg1
		var omg12 float64
		if (caps & capabilities.LongUnroll) != 0 {
			omg12 = E * (sig12 - (math.Atan2(ssig2, csig2) - math.Atan2(l.ssig1, l.csig1)) + (math.Atan2(E*somg2, comg2) - math.Atan2(E*l.somg1, l.comg1)))
		} else {
			omg12 = math.Atan2(somg2*l.comg1-comg2*l.somg1, comg2*l.comg1+somg2*l.somg1)
		}
		lam12 := omg12 + l.a3c*(sig12+(sinCosSeries(true, ssig2, csig2, l.c3a)-l.b31))
		lon12 := rad2deg(lam12)
		if (caps & capabilities.LongUnroll) != 0 {
			r.lon2 = l.lon1 + lon12
		} else {
			r.lon2 = angNormalize(angNormalize(l.lon1) + angNormalize(lon12))
		}
	}

	if (caps & capabilities.Latitude) != 0 {
		r.lat2 = atan2d(sbet2, l.g.f1*cbet2)
	}

	if (caps & capabilities.Azimuth) != 0 {
		r.azi2 = atan2d(salp2, calp2)
	}

	if (caps & (capabilities.ReducedLength | capabilities.GeodesicScale)) != 0 {
		b22 := sinCosSeries(true, ssig2, csig2, l.c2a)
		ab2 := (1 + l.a2m1) * (b22 - l.b21)
		j12 := (l.a1m1-l.a2m1)*sig12 + (ab1 - ab2)
		if (caps & capabilities.ReducedLength) != 0 {
			// Add parens around (l.csig1 * ssig2) and (l.ssig1 * csig2) to ensure
			// accurate cancellation in the case of coincident points.
			r.m12 = l.g.b * ((dn2*(l.csig1*ssig2) - l.dn1*(l.ssig1*csig2)) - l.csig1*csig2*j12)
		}
		if (caps & capabilities.GeodesicScale) != 0 {
			t := l.k2 * (ssig2 - l.ssig1) * (ssig2 + l.ssig1) / (l.dn1 + dn2)
			r.M12 = csig12 + (t*ssig2-csig2*j12)*l.ssig1/l.dn1
			r.M21 = csig12 - (t*l.ssig1-l.csig1*j12)*ssig2/dn2
		}
	}

	if (caps & capabilities.Area) != 0 {
		b42 := sinCosSeries(false, ssig2, csig2, l.c4a)
		var salp12, calp12 float64
		if l.calp0 == 0 || l.salp0 == 0 {
			// alp12 = alp2 - alp1, used in atan2 so no need to normalize
			salp12 = salp2*l.calp1 - calp2*l.salp1
			calp12 = calp2*l.calp1 + salp2*l.salp1
		} else {
			// tan(alp) = tan(alp0) * sec(sig)
			// tan(alp2-alp1) = (tan(alp2) -tan(alp1)) / (tan(alp2)*tan(alp1)+1)
			// = calp0 * salp0 * (csig1-csig2) / (salp0^2 + calp0^2 * csig1*csig2)
			// If csig12 > 0, write
			//   csig1 - csig2 = ssig12 * (csig1 * ssig12 / (1 + csig12) + ssig1)
			// else
			//   csig1 - csig2 = csig1 * (1 - csig12) + ssig12 * ssig1
			// No need to normalize
			var t float64
			if csig12 <= 0 {
				t = l.csig1*(1-csig12) + ssig12*l.ssig1
			} else {
				t = ssig12 * (l.csig1*ssig12/(1+csig12) + l.ssig1)
			}
			salp12 = l.calp0 * l.salp0 * t
			calp12 = sq(l.salp0) + sq(l.calp0)*l.csig1*csig2
		}
		r.S12 = l.g.c2*math.Atan2(salp12, calp12) + l.a4*(b42-l.b41)
	}
	return r
}

// SetDistance specifies the position of point 3 on the geodesic in terms of distance (meters). This
// is only useful if the Line instance was created with capabilities.DistanceIn.
func (l *Line) SetDistance(s13 float64) {
	l.s13 = s13
	l.a13 = l.PositionWithCapabilities(s13, capabilities.None).A12
}

// SetArc specifies the position of point 3 on the geodesic in terms of arc length (degrees). This
// is only useful if the Line instance was created with capabilities.Distance.
func (l *Line) SetArc(a13 float64) {
	l.a13 = a13
	l.s13 = l.ArcPositionWithCapabilities(a13, capabilities.Distance).S12
}

// Distance returns s13, the distance to point 3 (meters).
func (l *Line) Distance() float64 {
	return l.s13
}

// Arc returns a13, the arc length to point 3 (degrees).
func (l *Line) Arc() float64 {
	return l.a13
}
