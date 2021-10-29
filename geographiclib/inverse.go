package geographiclib

type inverseSolver struct {
	Geodesic
}

func (is *inverseSolver) inverse(lat1, lon1, lat2, lon2 float64, outmask int) *geodesicDataImpl {
	result := newGeodesicDataImpl()
	// Compute longitude difference (AngDiff does this carefully). Result is in [-180, 180] but -180
	// is only for west-going geodesics. 180 is for east-going and meridional geodesics.
	lat1, lat2 = latFix(lat1), latFix(lat2)
	result.lat1, result.lat2 = lat1, lat2
	// If really close to the equator, treat as on equator.
	lat1, lat2 = angRound(lat1), angRound(lat2)
	lon12, lon12s := angDiff(lon1, lon2)
	if (outmask & LongUnroll) != 0 {
		result.lon1 = lon1
		result.lon2 = (lon1 + lon12) + lon12s
	} else {
		result.lon1, result.lon2 = angNormalize(lon1), angNormalize(lon2)
	}
	// Make longitude difference positive.
	lonSign := sign(lon12)
	// If very close to being on the same half-meridian, then make it so.
	lon12 = lonSign * angRound(lon12)
	lon12s = angRound((180 - lon12) - lonSign*lon12s)
	lam12 := deg2rad(lon12)
	t := lon12
	if lon12 > 90 {
		t = lon12s
	}
	slam12, clam12 := sincosd(t)
	if lon12 > 90 {
		clam12 *= -1
	}

	// TODO: continue implementing me!

	lam12, slam12 = slam12, lam12 // this is nonsense
	return result
}
