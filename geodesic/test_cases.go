package geodesic

import (
	"fmt"
	"math"
	"testing"

	"github.com/pymaxion/geographiclib-go/geodesic/capabilities"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

var commonTestCases = []Data{
	{35.60777, -139.44815, 111.098748429560326,
		-11.17491, -69.95921, 129.289270889708762,
		8935244.5604818305, 80.50729714281974, 6273170.2055303837,
		0.16606318447386067, 0.16479116945612937, 12841384694976.432},
	{55.52454, 106.05087, 22.020059880982801,
		77.03196, 197.18234, 109.112041110671519,
		4105086.1713924406, 36.892740690445894, 3828869.3344387607,
		0.80076349608092607, 0.80101006984201008, 61674961290615.615},
	{-21.97856, 142.59065, -32.44456876433189,
		41.84138, 98.56635, -41.84359951440466,
		8394328.894657671, 75.62930491011522, 6161154.5773110616,
		0.24816339233950381, 0.24930251203627892, -6637997720646.717},
	{-66.99028, 112.2363, 173.73491240878403,
		-12.70631, 285.90344, 2.512956620913668,
		11150344.2312080241, 100.278634181155759, 6289939.5670446687,
		-0.17199490274700385, -0.17722569526345708, -121287239862139.744},
	{-17.42761, 173.34268, -159.033557661192928,
		-15.84784, 5.93557, -20.787484651536988,
		16076603.1631180673, 144.640108810286253, 3732902.1583877189,
		-0.81273638700070476, -0.81299800519154474, 97825992354058.708},
	{32.84994, 48.28919, 150.492927788121982,
		-56.28556, 202.29132, 48.113449399816759,
		16727068.9438164461, 150.565799985466607, 3147838.1910180939,
		-0.87334918086923126, -0.86505036767110637, -72445258525585.010},
	{6.96833, 52.74123, 92.581585386317712,
		-7.39675, 206.17291, 90.721692165923907,
		17102477.2496958388, 154.147366239113561, 2772035.6169917581,
		-0.89991282520302447, -0.89986892177110739, -1311796973197.995},
	{-50.56724, -16.30485, -105.439679907590164,
		-33.56571, -94.97412, -47.348547835650331,
		6455670.5118668696, 58.083719495371259, 5409150.7979815838,
		0.53053508035997263, 0.52988722644436602, 41071447902810.047},
	{-58.93002, -8.90775, 140.965397902500679,
		-8.91104, 133.13503, 19.255429433416599,
		11756066.0219864627, 105.755691241406877, 6151101.2270708536,
		-0.26548622269867183, -0.27068483874510741, -86143460552774.735},
	{-68.82867, -74.28391, 93.774347763114881,
		-50.63005, -8.36685, 34.65564085411343,
		3956936.926063544, 35.572254987389284, 3708890.9544062657,
		0.81443963736383502, 0.81420859815358342, -41845309450093.787},
	{-10.62672, -32.0898, -86.426713286747751,
		5.883, -134.31681, -80.473780971034875,
		11470869.3864563009, 103.387395634504061, 6184411.6622659713,
		-0.23138683500430237, -0.23155097622286792, 4198803992123.548},
	{-21.76221, 166.90563, 29.319421206936428,
		48.72884, 213.97627, 43.508671946410168,
		9098627.3986554915, 81.963476716121964, 6299240.9166992283,
		0.13965943368590333, 0.14152969707656796, 10024709850277.476},
	{-19.79938, -174.47484, 71.167275780171533,
		-11.99349, -154.35109, 65.589099775199228,
		2319004.8601169389, 20.896611684802389, 2267960.8703918325,
		0.93427001867125849, 0.93424887135032789, -3935477535005.785},
	{-11.95887, -116.94513, 92.712619830452549,
		4.57352, 7.16501, 78.64960934409585,
		13834722.5801401374, 124.688684161089762, 5228093.177931598,
		-0.56879356755666463, -0.56918731952397221, -9919582785894.853},
	{-87.85331, 85.66836, -65.120313040242748,
		66.48646, 16.09921, -4.888658719272296,
		17286615.3147144645, 155.58592449699137, 2635887.4729110181,
		-0.90697975771398578, -0.91095608883042767, 42667211366919.534},
	{1.74708, 128.32011, -101.584843631173858,
		-11.16617, 11.87109, -86.325793296437476,
		12942901.1241347408, 116.650512484301857, 5682744.8413270572,
		-0.44857868222697644, -0.44824490340007729, 10763055294345.653},
	{-25.72959, -144.90758, -153.647468693117198,
		-57.70581, -269.17879, -48.343983158876487,
		9413446.7452453107, 84.664533838404295, 6356176.6898881281,
		0.09492245755254703, 0.09737058264766572, 74515122850712.444},
	{-41.22777, 122.32875, 14.285113402275739,
		-7.57291, 130.37946, 10.805303085187369,
		3812686.035106021, 34.34330804743883, 3588703.8812128856,
		0.82605222593217889, 0.82572158200920196, -2456961531057.857},
	{11.01307, 138.25278, 79.43682622782374,
		6.62726, 247.05981, 103.708090215522657,
		11911190.819018408, 107.341669954114577, 6070904.722786735,
		-0.29767608923657404, -0.29785143390252321, 17121631423099.696},
	{-29.47124, 95.14681, -163.779130441688382,
		-27.46601, -69.15955, -15.909335945554969,
		13487015.8381145492, 121.294026715742277, 5481428.9945736388,
		-0.51527225545373252, -0.51556587964721788, 104679964020340.318},
}

type geodSolve struct {
	testNum     int
	description string
	logic       func(t *testing.T)
}

func (g *geodSolve) String() string {
	if g.description == "" {
		return fmt.Sprintf("GeodSolve%d", g.testNum)
	}
	return fmt.Sprintf("GeodSolve%d   %s", g.testNum, g.description)
}

var (
	geodSolve0 = geodSolve{
		testNum:     0,
		description: "",
		logic: func(t *testing.T) {
			r := WGS84.Inverse(40.6, -73.8, 49.01666667, 2.55)
			assert.InDelta(t, 53.47022, r.Azi1, 0.5e-5)
			assert.InDelta(t, 111.59367, r.Azi2, 0.5e-5)
			assert.InDelta(t, 5853226, r.S12, 0.5)
		},
	}

	geodSolve1 = geodSolve{
		testNum:     1,
		description: "",
		logic: func(t *testing.T) {
			r := WGS84.Direct(40.63972222, -73.77888889, 53.5, 5850e3)
			assert.InDelta(t, 49.01467, r.Lat2, 0.5e-5)
			assert.InDelta(t, 2.56106, r.Lon2, 0.5e-5)
			assert.InDelta(t, 111.62947, r.Azi2, 0.5e-5)
		},
	}

	geodSolve2 = geodSolve{
		testNum:     2,
		description: "Check fix for antipodal prolate bug found 2010-09-04",
		logic: func(t *testing.T) {
			geod, err := NewGeodesic(6.4e6, -1/150.0)
			require.Nil(t, err)

			r := geod.Inverse(0.07476, 0, -0.07476, 180)
			assert.InDelta(t, 90.00078, r.Azi1, 0.5e-5)
			assert.InDelta(t, 90.00078, r.Azi2, 0.5e-5)
			assert.InDelta(t, 20106193, r.S12, 0.5)

			r = geod.Inverse(0.1, 0, -0.1, 180)
			assert.InDelta(t, 90.00105, r.Azi1, 0.5e-5)
			assert.InDelta(t, 90.00105, r.Azi2, 0.5e-5)
			assert.InDelta(t, 20106193, r.S12, 0.5)
		},
	}

	geodSolve4 = geodSolve{
		testNum:     4,
		description: "Check fix for short line bug found 2010-05-21",
		logic: func(t *testing.T) {
			r := WGS84.Inverse(36.493349428792, 0, 36.49334942879201, .0000008)
			assert.InDelta(t, 0.072, r.S12, 0.5e-3)
		},
	}

	geodSolve5 = geodSolve{
		testNum:     5,
		description: "Check fix for point2=pole bug found 2010-05-03",
		logic: func(t *testing.T) {
			r := WGS84.Direct(0.01777745589997, 30, 0, 10e6)
			assert.InDelta(t, 90, r.Lat2, 0.5e-5)
			if r.Lon2 < 0 {
				assert.InDelta(t, -150, r.Lon2, 0.5e-5)
				assert.InDelta(t, 180, math.Abs(r.Azi2), 0.5e-5)
			} else {
				assert.InDelta(t, 30, r.Lon2, 0.5e-5)
				assert.InDelta(t, 0, r.Azi2, 0.5e-5)
			}
		},
	}

	geodSolve6 = geodSolve{
		testNum:     6,
		description: "Check fix for volatile sbet12a bug found 2011-06-25 (gcc 4.4.4 x86-O3). Found again on 2012-03-27 with tdm-mingw32 (g++ 4.6.1).",
		logic: func(t *testing.T) {
			r := WGS84.Inverse(88.202499451857, 0, -88.202499451857, 179.981022032992859592)
			assert.InDelta(t, 20003898.214, r.S12, 0.5e-3)

			r = WGS84.Inverse(89.262080389218, 0, -89.262080389218, 179.992207982775375662)
			assert.InDelta(t, 20003925.854, r.S12, 0.5e-3)

			r = WGS84.Inverse(89.333123580033, 0, -89.333123580032997687, 179.99295812360148422)
			assert.InDelta(t, 20003926.881, r.S12, 0.5e-3)
		},
	}

	geodSolve9 = geodSolve{
		testNum:     9,
		description: "Check fix for volatile x bug found 2011-06-25 (gcc 4.4.4 x86 -O3)",
		logic: func(t *testing.T) {
			r := WGS84.Inverse(56.320923501171, 0, -56.320923501171, 179.664747671772880215)
			assert.InDelta(t, 19993558.287, r.S12, 0.5e-3)
		},
	}

	geodSolve10 = geodSolve{
		testNum:     10,
		description: "Check fix for adjust tol1_ bug found 2011-06-25 (Visual Studio 10 rel + debug)",
		logic: func(t *testing.T) {
			r := WGS84.Inverse(52.784459512564, 0, -52.784459512563990912, 179.634407464943777557)
			assert.InDelta(t, 19991596.095, r.S12, 0.5e-3)
		},
	}

	geodSolve11 = geodSolve{
		testNum:     11,
		description: "Check fix for bet2 = -bet1 bug found 2011-06-25 (Visual Studio 10 rel + debug)",
		logic: func(t *testing.T) {
			r := WGS84.Inverse(48.522876735459, 0, -48.52287673545898293, 179.599720456223079643)
			assert.InDelta(t, 19989144.774, r.S12, 0.5e-3)
		},
	}

	geodSolve12 = geodSolve{
		testNum:     12,
		description: "Check fix for inverse geodesics on extreme prolate/oblate ellipsoids. Reported 2012-08-29 Stefan Guenther <stefan.gunther@embl.de>; fixed 2012-10-07",
		logic: func(t *testing.T) {
			geod, err := NewGeodesic(89.8, -1.83)
			require.Nil(t, err)

			r := geod.Inverse(0, 0, -10, 160)
			assert.InDelta(t, 120.27, r.Azi1, 1e-2)
			assert.InDelta(t, 105.15, r.Azi2, 1e-2)
			assert.InDelta(t, 266.7, r.S12, 1e-1)
		},
	}

	geodSolve14 = geodSolve{
		testNum:     14,
		description: "Check fix for inverse ignoring lon12 = nan",
		logic: func(t *testing.T) {
			r := WGS84.Inverse(0, 0, 1, math.NaN())
			assert.True(t, math.IsNaN(r.Azi1))
			assert.True(t, math.IsNaN(r.Azi2))
			assert.True(t, math.IsNaN(r.S12))
		},
	}

	geodSolve15 = geodSolve{
		testNum:     15,
		description: "Initial implementation of Math::eatanhe was wrong for e^2 < 0. This checks that this is fixed.",
		logic: func(t *testing.T) {
			geod, err := NewGeodesic(6.4e6, -1/150.0)
			require.Nil(t, err)

			r := geod.DirectWithCapabilities(1, 2, 3, 4, capabilities.Area)
			assert.InDelta(t, 23700, r.S12Area, 0.5)
		},
	}

	geodSolve17 = geodSolve{
		testNum:     17,
		description: "Check fix for LONG_UNROLL bug found on 2015-05-07",
		logic: func(t *testing.T) {
			r := WGS84.DirectWithCapabilities(40, -75, -10, 2e7, capabilities.Standard|capabilities.LongUnroll)
			assert.InDelta(t, -39, r.Lat2, 1)
			assert.InDelta(t, -254, r.Lon2, 1)
			assert.InDelta(t, -170, r.Azi2, 1)

			line := WGS84.Line(40, -75, -10)
			r = line.PositionWithCapabilities(2e7, capabilities.Standard|capabilities.LongUnroll)
			assert.InDelta(t, -39, r.Lat2, 1)
			assert.InDelta(t, -254, r.Lon2, 1)
			assert.InDelta(t, -170, r.Azi2, 1)

			r = WGS84.Direct(40, -75, -10, 2e7)
			assert.InDelta(t, -39, r.Lat2, 1)
			assert.InDelta(t, 105, r.Lon2, 1)
			assert.InDelta(t, -170, r.Azi2, 1)

			r = line.Position(2e7)
			assert.InDelta(t, -39, r.Lat2, 1)
			assert.InDelta(t, 105, r.Lon2, 1)
			assert.InDelta(t, -170, r.Azi2, 1)
		},
	}

	geodSolve26 = geodSolve{
		testNum:     26,
		description: "Check 0/0 problem with area calculation on sphere 2015-09-08",
		logic: func(t *testing.T) {
			geod, err := NewGeodesic(6.4e6, 0)
			require.Nil(t, err)

			r := geod.InverseWithCapabilities(1, 2, 3, 4, capabilities.Area)
			assert.InDelta(t, 49911046115.0, r.S12Area, 0.5)
		},
	}

	geodSolve28 = geodSolve{
		testNum:     28,
		description: "Check for bad placement of assignment of r.a12 with |f| > 0.01 (bug in Java implementation fixed on 2015-05-19).",
		logic: func(t *testing.T) {
			geod, err := NewGeodesic(6.4e6, 0.1)
			require.Nil(t, err)

			r := geod.Direct(1, 2, 10, 5e6)
			assert.InDelta(t, 48.55570690, r.A12, 0.5e-8)
		},
	}

	geodSolve29 = geodSolve{
		testNum:     29,
		description: "Check longitude unrolling with inverse calculation 2015-09-16",
		logic: func(t *testing.T) {
			r := WGS84.Inverse(0, 539, 0, 181)
			assert.InDelta(t, 179, r.Lon1, 1e-10)
			assert.InDelta(t, -179, r.Lon2, 1e-10)
			assert.InDelta(t, 222639, r.S12, 0.5)

			r = WGS84.InverseWithCapabilities(0, 539, 0, 181, capabilities.Standard|capabilities.LongUnroll)
			assert.InDelta(t, 539, r.Lon1, 1e-10)
			assert.InDelta(t, 541, r.Lon2, 1e-10)
			assert.InDelta(t, 222639, r.S12, 0.5)
		},
	}

	geodSolve33 = geodSolve{
		testNum:     33,
		description: "Check max(-0.0,+0.0) issues 2015-08-22 (triggered by bugs in Octave -- sind(-0.0) = +0.0 -- and in some version of Visual Studio -- fmod(-0.0, 360.0) = +0.0)",
		logic: func(t *testing.T) {
			r := WGS84.Inverse(0, 0, 0, 179)
			assert.InDelta(t, 90.00000, r.Azi1, 0.5e-5)
			assert.InDelta(t, 90.00000, r.Azi2, 0.5e-5)
			assert.InDelta(t, 19926189, r.S12, 0.5)

			r = WGS84.Inverse(0, 0, 0, 179.5)
			assert.InDelta(t, 55.96650, r.Azi1, 0.5e-5)
			assert.InDelta(t, 124.03350, r.Azi2, 0.5e-5)
			assert.InDelta(t, 19980862, r.S12, 0.5)

			r = WGS84.Inverse(0, 0, 0, 180)
			assert.InDelta(t, 0.00000, r.Azi1, 0.5e-5)
			assert.InDelta(t, 180.00000, math.Abs(r.Azi2), 0.5e-5)
			assert.InDelta(t, 20003931, r.S12, 0.5)

			r = WGS84.Inverse(0, 0, 1, 180)
			assert.InDelta(t, 0.00000, r.Azi1, 0.5e-5)
			assert.InDelta(t, 180.00000, math.Abs(r.Azi2), 0.5e-5)
			assert.InDelta(t, 19893357, r.S12, 0.5)

			geod, err := NewGeodesic(6.4e6, 0)
			require.Nil(t, err)

			r = geod.Inverse(0, 0, 0, 179)
			assert.InDelta(t, 90.00000, r.Azi1, 0.5e-5)
			assert.InDelta(t, 90.00000, r.Azi2, 0.5e-5)
			assert.InDelta(t, 19994492, r.S12, 0.5)

			r = geod.Inverse(0, 0, 0, 180)
			assert.InDelta(t, 0.00000, r.Azi1, 0.5e-5)
			assert.InDelta(t, 180.00000, math.Abs(r.Azi2), 0.5e-5)
			assert.InDelta(t, 20106193, r.S12, 0.5)

			r = geod.Inverse(0, 0, 1, 180)
			assert.InDelta(t, 0.00000, r.Azi1, 0.5e-5)
			assert.InDelta(t, 180.00000, math.Abs(r.Azi2), 0.5e-5)
			assert.InDelta(t, 19994492, r.S12, 0.5)

			geod, err = NewGeodesic(6.4e6, -1/300.0)
			require.Nil(t, err)

			r = geod.Inverse(0, 0, 0, 179)
			assert.InDelta(t, 90.00000, r.Azi1, 0.5e-5)
			assert.InDelta(t, 90.00000, r.Azi2, 0.5e-5)
			assert.InDelta(t, 19994492, r.S12, 0.5)

			r = geod.Inverse(0, 0, 0, 180)
			assert.InDelta(t, 90.00000, r.Azi1, 0.5e-5)
			assert.InDelta(t, 90.00000, r.Azi2, 0.5e-5)
			assert.InDelta(t, 20106193, r.S12, 0.5)

			r = geod.Inverse(0, 0, 0.5, 180)
			assert.InDelta(t, 33.02493, r.Azi1, 0.5e-5)
			assert.InDelta(t, 146.97364, r.Azi2, 0.5e-5)
			assert.InDelta(t, 20082617, r.S12, 0.5)

			r = geod.Inverse(0, 0, 1, 180)
			assert.InDelta(t, 0.00000, r.Azi1, 0.5e-5)
			assert.InDelta(t, 180.00000, math.Abs(r.Azi2), 0.5e-5)
			assert.InDelta(t, 20027270, r.S12, 0.5)
		},
	}

	geodSolve55 = geodSolve{
		testNum:     55,
		description: "Check fix for nan + point on equator or pole not returning all nans in Geodesic::Inverse, found 2015-09-23.",
		logic: func(t *testing.T) {
			r := WGS84.Inverse(math.NaN(), 0, 0, 90)
			assert.True(t, math.IsNaN(r.Azi1))
			assert.True(t, math.IsNaN(r.Azi2))
			assert.True(t, math.IsNaN(r.S12))

			r = WGS84.Inverse(math.NaN(), 0, 90, 3)
			assert.True(t, math.IsNaN(r.Azi1))
			assert.True(t, math.IsNaN(r.Azi2))
			assert.True(t, math.IsNaN(r.S12))
		},
	}

	geodSolve59 = geodSolve{
		testNum:     59,
		description: "Check for points close with longitudes close to 180 deg apart.",
		logic: func(t *testing.T) {
			r := WGS84.Inverse(5, 0.00000000000001, 10, 180)
			assert.InDelta(t, 0.000000000000035, r.Azi1, 1.5e-14)
			assert.InDelta(t, 179.99999999999996, r.Azi2, 1.5e-14)
			assert.InDelta(t, 18345191.174332713, r.S12, 5e-9)
		},
	}

	geodSolve61 = geodSolve{
		testNum:     61,
		description: "Make sure small negative azimuths are west-going",
		logic: func(t *testing.T) {
			r := WGS84.DirectWithCapabilities(45, 0, -0.000000000000000003, 1e7, capabilities.Standard|capabilities.LongUnroll)
			assert.InDelta(t, 45.30632, r.Lat2, 0.5e-5)
			assert.InDelta(t, -180, r.Lon2, 0.5e-5)
			assert.InDelta(t, 180, math.Abs(r.Azi2), 0.5e-5)

			l := WGS84.InverseLine(45, 0, 80, -0.000000000000000003)
			r = l.PositionWithCapabilities(1e7, capabilities.Standard|capabilities.LongUnroll)
			assert.InDelta(t, 45.30632, r.Lat2, 0.5e-5)
			assert.InDelta(t, -180, r.Lon2, 0.5e-5)
			assert.InDelta(t, 180, math.Abs(r.Azi2), 0.5e-5)
		},
	}

	geodSolve65 = geodSolve{
		testNum:     65,
		description: "Check for bug in east-going check in GeodesicLine (needed to check for sign of 0) and sign error in area calculation due to a bogus override of the code for alp12.  Found/fixed on 2015-12-19.",
		logic: func(t *testing.T) {
			l := WGS84.InverseLine(30, -0.000000000000000001, -31, 180)
			r := l.PositionWithCapabilities(1e7, capabilities.All|capabilities.LongUnroll)
			assert.InDelta(t, 30.00000, r.Lat1, 0.5e-5)
			assert.InDelta(t, -0.00000, r.Lon1, 0.5e-5)
			assert.InDelta(t, 180.00000, math.Abs(r.Azi1), 0.5e-5)
			assert.InDelta(t, -60.23169, r.Lat2, 0.5e-5)
			assert.InDelta(t, -0.00000, r.Lon2, 0.5e-5)
			assert.InDelta(t, 180.00000, math.Abs(r.Azi2), 0.5e-5)
			assert.InDelta(t, 10000000, r.S12, 0.5)
			assert.InDelta(t, 90.06544, r.A12, 0.5e-5)
			assert.InDelta(t, 6363636, r.M12Reduced, 0.5)
			assert.InDelta(t, -0.0012834, r.M12, 0.5e7)
			assert.InDelta(t, 0.0013749, r.M21, 0.5e-7)
			assert.InDelta(t, 0, r.S12Area, 0.5)

			r = l.PositionWithCapabilities(2e7, capabilities.All|capabilities.LongUnroll)
			assert.InDelta(t, 30.00000, r.Lat1, 0.5e-5)
			assert.InDelta(t, -0.00000, r.Lon1, 0.5e-5)
			assert.InDelta(t, 180.00000, math.Abs(r.Azi1), 0.5e-5)
			assert.InDelta(t, -30.03547, r.Lat2, 0.5e-5)
			assert.InDelta(t, -180.00000, r.Lon2, 0.5e-5)
			assert.InDelta(t, -0.00000, r.Azi2, 0.5e-5)
			assert.InDelta(t, 20000000, r.S12, 0.5)
			assert.InDelta(t, 179.96459, r.A12, 0.5e-5)
			assert.InDelta(t, 54342, r.M12Reduced, 0.5)
			assert.InDelta(t, -1.0045592, r.M12, 0.5e7)
			assert.InDelta(t, -0.9954339, r.M21, 0.5e-7)
			assert.InDelta(t, 127516405431022.0, r.S12Area, 0.5)
		},
	}

	geodSolve69 = geodSolve{
		testNum:     69,
		description: "Check for InverseLine if line is slightly west of S and that s13 is correctly set.",
		logic: func(t *testing.T) {
			l := WGS84.InverseLine(-5, -0.000000000000002, -10, 180)
			r := l.PositionWithCapabilities(2e7, capabilities.Standard|capabilities.LongUnroll)
			assert.InDelta(t, 4.96445, r.Lat2, 0.5e-5)
			assert.InDelta(t, -180.00000, r.Lon2, 0.5e-5)
			assert.InDelta(t, -0.00000, r.Azi2, 0.5e-5)

			r = l.PositionWithCapabilities(0.5*l.Distance(), capabilities.Standard|capabilities.LongUnroll)
			assert.InDelta(t, -87.52461, r.Lat2, 0.5e-5)
			assert.InDelta(t, -0.00000, r.Lon2, 0.5e-5)
			assert.InDelta(t, -180.00000, r.Azi2, 0.5e-5)
		},
	}

	geodSolve71 = geodSolve{
		testNum:     71,
		description: "Check that DirectLine sets s13.",
		logic: func(t *testing.T) {
			l := WGS84.DirectLine(1, 2, 45, 1e7)
			r := l.PositionWithCapabilities(0.5*l.Distance(), capabilities.Standard|capabilities.LongUnroll)
			assert.InDelta(t, 30.92625, r.Lat2, 0.5e-5)
			assert.InDelta(t, 37.54640, r.Lon2, 0.5e-5)
			assert.InDelta(t, 55.43104, r.Azi2, 0.5e-5)
		},
	}

	geodSolve73 = geodSolve{
		testNum:     73,
		description: "Check for backwards from the pole bug reported by Anon on 2016-02-13. This only affected the Java implementation. It was introduced in Java version 1.44 and fixed in 1.46-SNAPSHOT on 2016-01-17. Also the + sign on azi2 is a check on the normalizing of azimuths (converting -0.0 to +0.0).",
		logic: func(t *testing.T) {
			r := WGS84.Direct(90, 10, 180, -1e6)
			assert.InDelta(t, 81.04623, r.Lat2, 0.5e-5)
			assert.InDelta(t, -170, r.Lon2, 0.5e-5)
			assert.InDelta(t, 0, r.Azi2, 0.5e-5)
			assert.False(t, math.Signbit(r.Azi2))
		},
	}

	geodSolve74 = geodSolve{
		testNum:     74,
		description: "Check fix for inaccurate areas, bug introduced in v1.46, fixed 2015-10-16.",
		logic: func(t *testing.T) {
			r := WGS84.InverseWithCapabilities(54.1589, 15.3872, 54.1591, 15.3877, capabilities.All)
			assert.InDelta(t, 55.723110355, r.Azi1, 5e-9)
			assert.InDelta(t, 55.723515675, r.Azi2, 5e-9)
			assert.InDelta(t, 39.527686385, r.S12, 5e-9)
			assert.InDelta(t, 0.000355495, r.A12, 5e-9)
			assert.InDelta(t, 39.527686385, r.M12Reduced, 5e-9)
			assert.InDelta(t, 0.999999995, r.M12, 5e-9)
			assert.InDelta(t, 0.999999995, r.M21, 5e-9)
			assert.InDelta(t, 286698586.30197, r.S12Area, 5e-4)
		},
	}

	geodSolve76 = geodSolve{
		testNum:     76,
		description: "The distance from Wellington and Salamanca (a classic failure of Vincenty)",
		logic: func(t *testing.T) {
			r := WGS84.Inverse(-(41 + 19/60.0), 174+49/60.0, 40+58/60.0, -(5 + 30/60.0))
			assert.InDelta(t, 160.39137649664, r.Azi1, 0.5e-11)
			assert.InDelta(t, 19.50042925176, r.Azi2, 0.5e-11)
			assert.InDelta(t, 19960543.857179, r.S12, 0.5e-6)
		},
	}

	geodSolve78 = geodSolve{
		testNum:     78,
		description: "An example where the NGS calculator fails to converge",
		logic: func(t *testing.T) {
			r := WGS84.Inverse(27.2, 0.0, -27.1, 179.5)
			assert.InDelta(t, 45.82468716758, r.Azi1, 0.5e-11)
			assert.InDelta(t, 134.22776532670, r.Azi2, 0.5e-11)
			assert.InDelta(t, 19974354.765767, r.S12, 0.5e-6)

		},
	}

	geodSolve80 = geodSolve{
		testNum:     80,
		description: "Some tests to add code coverage: computing scale in special cases + zero length geodesic (includes GeodSolve80 - GeodSolve83).",
		logic: func(t *testing.T) {
			r := WGS84.InverseWithCapabilities(0, 0, 0, 90, capabilities.GeodesicScale)
			assert.InDelta(t, -0.00528427534, r.M12, 0.5e-10)
			assert.InDelta(t, -0.00528427534, r.M21, 0.5e-10)

			r = WGS84.InverseWithCapabilities(0, 0, 1e-6, 1e-6, capabilities.GeodesicScale)
			assert.InDelta(t, 1, r.M12, 0.5e-10)
			assert.InDelta(t, 1, r.M21, 0.5e-10)

			r = WGS84.InverseWithCapabilities(20.001, 0, 20.001, 0, capabilities.All)
			assert.InDelta(t, 0, r.A12, 1e-13)
			assert.InDelta(t, 0, r.S12, 1e-8)
			assert.InDelta(t, 180, r.Azi1, 1e-13)
			assert.InDelta(t, 180, r.Azi2, 1e-13)
			assert.InDelta(t, 0, r.M12Reduced, 1e-8)
			assert.InDelta(t, 1, r.M12, 1e-15)
			assert.InDelta(t, 1, r.M21, 1e-15)
			assert.InDelta(t, 0, r.S12Area, 1e-10)
			assert.False(t, math.Signbit(r.A12))
			assert.False(t, math.Signbit(r.S12))
			assert.False(t, math.Signbit(r.M12Reduced))

			r = WGS84.InverseWithCapabilities(90, 0, 90, 180, capabilities.All)
			assert.InDelta(t, 0, r.A12, 1e-13)
			assert.InDelta(t, 0, r.S12, 1e-8)
			assert.InDelta(t, 0, r.Azi1, 1e-13)
			assert.InDelta(t, 180, r.Azi2, 1e-13)
			assert.InDelta(t, 0, r.M12Reduced, 1e-8)
			assert.InDelta(t, 1, r.M12, 1e-15)
			assert.InDelta(t, 1, r.M21, 1e-15)
			assert.InDelta(t, 127516405431022.0, r.S12Area, 0.5)

			// An incapable line which can't take distance as input
			line := WGS84.LineWithCapabilities(1, 2, 90, capabilities.Latitude)
			r = line.PositionWithCapabilities(1000, capabilities.None)
			assert.True(t, math.IsNaN(r.A12))
		},
	}

	geodSolve84 = geodSolve{
		testNum:     84,
		description: "Tests for python implementation to check fix for range errors with {fmod,sin,cos}(inf) (includes GeodSolve84 - GeodSolve91).",
		logic: func(t *testing.T) {
			r := WGS84.Direct(0, 0, 90, math.Inf(1))
			assert.True(t, math.IsNaN(r.Lat2))
			assert.True(t, math.IsNaN(r.Lon2))
			assert.True(t, math.IsNaN(r.Azi2))

			r = WGS84.Direct(0, 0, 90, math.NaN())
			assert.True(t, math.IsNaN(r.Lat2))
			assert.True(t, math.IsNaN(r.Lon2))
			assert.True(t, math.IsNaN(r.Azi2))

			r = WGS84.Direct(0, 0, math.Inf(1), 1000)
			assert.True(t, math.IsNaN(r.Lat2))
			assert.True(t, math.IsNaN(r.Lon2))
			assert.True(t, math.IsNaN(r.Azi2))

			r = WGS84.Direct(0, 0, math.NaN(), 1000)
			assert.True(t, math.IsNaN(r.Lat2))
			assert.True(t, math.IsNaN(r.Lon2))
			assert.True(t, math.IsNaN(r.Azi2))

			r = WGS84.Direct(0, math.Inf(1), 90, 1000)
			assert.True(t, r.Lat2 == 0)
			assert.True(t, math.IsNaN(r.Lon2))
			assert.True(t, r.Azi2 == 90)

			r = WGS84.Direct(0, math.NaN(), 90, 1000)
			assert.True(t, r.Lat2 == 0)
			assert.True(t, math.IsNaN(r.Lon2))
			assert.True(t, r.Azi2 == 90)

			r = WGS84.Direct(math.Inf(1), 0, 90, 1000)
			assert.True(t, math.IsNaN(r.Lat2))
			assert.True(t, math.IsNaN(r.Lon2))
			assert.True(t, math.IsNaN(r.Azi2))

			r = WGS84.Direct(math.NaN(), 0, 90, 1000)
			assert.True(t, math.IsNaN(r.Lat2))
			assert.True(t, math.IsNaN(r.Lon2))
			assert.True(t, math.IsNaN(r.Azi2))
		},
	}

	geodSolve92 = geodSolve{
		testNum:     92,
		description: "Check fix for inaccurate hypot with python 3.[89]. Problem reported by agdhruv https://github.com/geopy/geopy/issues/466; see https://bugs.python.org/issue43088",
		logic: func(t *testing.T) {
			r := WGS84.Inverse(37.757540000000006, -122.47018, 37.75754, -122.470177)
			assert.InDelta(t, 89.99999923, r.Azi1, 1e-7)
			assert.InDelta(t, 90.00000106, r.Azi2, 1e-7)
			assert.InDelta(t, 0.264, r.S12, 0.5e-3)
		},
	}
)

type planimeterTest struct {
	testNum     int
	description string
	logic       func(t *testing.T)
}

func (p *planimeterTest) String() string {
	if p.description == "" {
		return fmt.Sprintf("Planimeter%d", p.testNum)
	}
	return fmt.Sprintf("Planimeter%d   %s", p.testNum, p.description)
}

func planimeter(points [][]float64) PolygonResult {
	polygon := NewPolygonArea(WGS84, false)
	for i := 0; i < len(points); i++ {
		polygon.AddPoint(points[i][0], points[i][1])
	}
	return polygon.Compute(false, true)
}

func polyLength(points [][]float64) PolygonResult {
	polyline := NewPolygonArea(WGS84, true)
	for i := 0; i < len(points); i++ {
		polyline.AddPoint(points[i][0], points[i][1])
	}
	return polyline.Compute(false, true)
}

var (
	planimeter0 = planimeterTest{
		testNum:     0,
		description: "Check fix for pole-encircling bug found 2011-03-16",
		logic: func(t *testing.T) {
			pa := [][]float64{{89, 0}, {89, 90}, {89, 180}, {89, 270}}
			a := planimeter(pa)
			assert.InDelta(t, 631819.8745, a.Perimeter, 1e-4)
			assert.InDelta(t, 24952305678.0, a.Area, 1)

			pb := [][]float64{{-89, 0}, {-89, 90}, {-89, 180}, {-89, 270}}
			a = planimeter(pb)
			assert.InDelta(t, 631819.8745, a.Perimeter, 1e-4)
			assert.InDelta(t, -24952305678.0, a.Area, 1)

			pc := [][]float64{{0, -1}, {-1, 0}, {0, 1}, {1, 0}}
			a = planimeter(pc)
			assert.InDelta(t, 627598.2731, a.Perimeter, 1e-4)
			assert.InDelta(t, 24619419146.0, a.Area, 1)

			pd := [][]float64{{90, 0}, {0, 0}, {0, 90}}
			a = planimeter(pd)
			assert.InDelta(t, 30022685, a.Perimeter, 1)
			assert.InDelta(t, 63758202715511.0, a.Area, 1)

			a = polyLength(pd)
			assert.InDelta(t, 20020719, a.Perimeter, 1)
			assert.True(t, math.IsNaN(a.Area))
		},
	}

	planimeter5 = planimeterTest{
		testNum:     5,
		description: "Check fix for planimeter pole crossing bug found 2011-06-24",
		logic: func(t *testing.T) {
			points := [][]float64{{89, 0.1}, {89, 90.1}, {89, -179.9}}
			a := planimeter(points)
			assert.InDelta(t, 539297, a.Perimeter, 1)
			assert.InDelta(t, 12476152838.5, a.Area, 1)
		},
	}

	planimeter6 = planimeterTest{
		testNum:     6,
		description: "Check fix for Planimeter lon12 rounding bug found 2012-12-03",
		logic: func(t *testing.T) {
			pa := [][]float64{{9, -0.00000000000001}, {9, 180}, {9, 0}}
			a := planimeter(pa)
			assert.InDelta(t, 36026861, a.Perimeter, 1)
			assert.InDelta(t, 0, a.Area, 1)

			pb := [][]float64{{9, 0.00000000000001}, {9, 0}, {9, 180}}
			a = planimeter(pb)
			assert.InDelta(t, 36026861, a.Perimeter, 1)
			assert.InDelta(t, 0, a.Area, 1)

			pc := [][]float64{{9, 0.00000000000001}, {9, 180}, {9, 0}}
			a = planimeter(pc)
			assert.InDelta(t, 36026861, a.Perimeter, 1)
			assert.InDelta(t, 0, a.Area, 1)

			pd := [][]float64{{9, -0.00000000000001}, {9, 0}, {9, 180}}
			a = planimeter(pd)
			assert.InDelta(t, 36026861, a.Perimeter, 1)
			assert.InDelta(t, 0, a.Area, 1)
		},
	}

	planimeter12 = planimeterTest{
		testNum:     12,
		description: "Area of arctic circle (not really -- adjunct to rhumb-area test)",
		logic: func(t *testing.T) {
			points := [][]float64{{66.562222222, 0}, {66.562222222, 180}}
			a := planimeter(points)
			assert.InDelta(t, 10465729, a.Perimeter, 1)
			assert.InDelta(t, 0, a.Area, 1)
		},
	}

	planimeter13 = planimeterTest{
		testNum:     13,
		description: "Check encircling pole twice",
		logic: func(t *testing.T) {
			points := [][]float64{{89, -360}, {89, -240}, {89, -120}, {89, 0}, {89, 120}, {89, 240}}
			a := planimeter(points)
			assert.InDelta(t, 1160741, a.Perimeter, 1)
			assert.InDelta(t, 32415230256.0, a.Area, 1)
		},
	}

	planimeter15 = planimeterTest{
		testNum:     15,
		description: "Coverage tests, includes Planimeter15 - Planimeter18 (combinations of reverse and sign) + calls to testpoint, testedge.",
		logic: func(t *testing.T) {
			lat := []float64{2, 1, 3}
			lon := []float64{1, 2, 3}
			r := 18454562325.45119
			a0 := 510065621724088.5093 // ellipsoid area

			polygon := NewPolygonArea(WGS84, false)
			polygon.AddPoint(lat[0], lon[0])
			polygon.AddPoint(lat[1], lon[1])

			a := polygon.TestPoint(lat[2], lon[2], false, true)
			assert.InDelta(t, r, a.Area, 0.5)
			a = polygon.TestPoint(lat[2], lon[2], false, false)
			assert.InDelta(t, r, a.Area, 0.5)
			a = polygon.TestPoint(lat[2], lon[2], true, true)
			assert.InDelta(t, -r, a.Area, 0.5)
			a = polygon.TestPoint(lat[2], lon[2], true, false)
			assert.InDelta(t, a0-r, a.Area, 0.5)

			inv := WGS84.Inverse(lat[1], lon[1], lat[2], lon[2])

			a = polygon.TestEdge(inv.Azi1, inv.S12, false, true)
			assert.InDelta(t, r, a.Area, 0.5)
			a = polygon.TestEdge(inv.Azi1, inv.S12, false, false)
			assert.InDelta(t, r, a.Area, 0.5)
			a = polygon.TestEdge(inv.Azi1, inv.S12, true, true)
			assert.InDelta(t, -r, a.Area, 0.5)
			a = polygon.TestEdge(inv.Azi1, inv.S12, true, false)
			assert.InDelta(t, a0-r, a.Area, 0.5)

			polygon.AddPoint(lat[2], lon[2])

			a = polygon.Compute(false, true)
			assert.InDelta(t, r, a.Area, 0.5)
			a = polygon.Compute(false, false)
			assert.InDelta(t, r, a.Area, 0.5)
			a = polygon.Compute(true, true)
			assert.InDelta(t, -r, a.Area, 0.5)
			a = polygon.Compute(true, false)
			assert.InDelta(t, a0-r, a.Area, 0.5)
		},
	}

	planimeter19 = planimeterTest{
		testNum:     19,
		description: "Coverage tests, includes Planimeter19 - Planimeter20 (degenerate polygons) + extra cases.",
		logic: func(t *testing.T) {
			polygon := NewPolygonArea(WGS84, false)
			a := polygon.Compute(false, true)
			assert.True(t, a.Area == 0)
			assert.True(t, a.Perimeter == 0)
			a = polygon.TestPoint(1, 1, false, true)
			assert.True(t, a.Area == 0)
			assert.True(t, a.Perimeter == 0)
			a = polygon.TestEdge(90, 1000, false, true)
			assert.True(t, math.IsNaN(a.Area))
			assert.True(t, math.IsNaN(a.Perimeter))
			polygon.AddPoint(1, 1)
			a = polygon.Compute(false, true)
			assert.True(t, a.Area == 0)
			assert.True(t, a.Perimeter == 0)

			polyline := NewPolygonArea(WGS84, true)
			a = polyline.Compute(false, true)
			assert.True(t, a.Perimeter == 0)
			a = polyline.TestPoint(1, 1, false, true)
			assert.True(t, a.Perimeter == 0)
			a = polyline.TestEdge(90, 1000, false, true)
			assert.True(t, math.IsNaN(a.Perimeter))
			polyline.AddPoint(1, 1)
			a = polyline.Compute(false, true)
			assert.True(t, a.Perimeter == 0)
			polygon.AddPoint(1, 1)
			a = polyline.TestEdge(90, 1000, false, true)
			assert.InDelta(t, 1000, a.Perimeter, 1e-10)
			a = polyline.TestPoint(2, 2, false, true)
			assert.InDelta(t, 156876.149, a.Perimeter, 0.5e-3)
		},
	}

	planimeter21 = planimeterTest{
		testNum:     21,
		description: "Some test to add code coverage: multiple circlings of pole (includes Planimeter21 - Planimeter28) + invocations via testpoint and testedge.",
		logic: func(t *testing.T) {
			lat := 45.
			azi := 39.2144607176828184218
			s := 8420705.40957178156285
			r := 39433884866571.4277   // Area for one circuit
			a0 := 510065621724088.5093 // Ellipsoid area
			polygon := NewPolygonArea(WGS84, false)
			polygon.AddPoint(lat, 60)
			polygon.AddPoint(lat, 180)
			polygon.AddPoint(lat, -60)
			polygon.AddPoint(lat, 60)
			polygon.AddPoint(lat, 180)
			polygon.AddPoint(lat, -60)

			for i := 3.; i <= 4.; i++ {
				polygon.AddPoint(lat, 60)
				polygon.AddPoint(lat, 180)
				a := polygon.TestPoint(lat, -60, false, true)
				assert.InDelta(t, i*r, a.Area, 0.5)
				a = polygon.TestPoint(lat, -60, false, false)
				assert.InDelta(t, i*r, a.Area, 0.5)
				a = polygon.TestPoint(lat, -60, true, true)
				assert.InDelta(t, -i*r, a.Area, 0.5)
				a = polygon.TestPoint(lat, -60, true, false)
				assert.InDelta(t, -i*r+a0, a.Area, 0.5)
				a = polygon.TestEdge(azi, s, false, true)
				assert.InDelta(t, i*r, a.Area, 0.5)
				a = polygon.TestEdge(azi, s, false, false)
				assert.InDelta(t, i*r, a.Area, 0.5)
				a = polygon.TestEdge(azi, s, true, true)
				assert.InDelta(t, -i*r, a.Area, 0.5)
				a = polygon.TestEdge(azi, s, true, false)
				assert.InDelta(t, -i*r+a0, a.Area, 0.5)
				polygon.AddPoint(lat, -60)
				a = polygon.Compute(false, true)
				assert.InDelta(t, i*r, a.Area, 0.5)
				a = polygon.Compute(false, false)
				assert.InDelta(t, i*r, a.Area, 0.5)
				a = polygon.Compute(true, true)
				assert.InDelta(t, -i*r, a.Area, 0.5)
				a = polygon.Compute(true, false)
				assert.InDelta(t, -i*r+a0, a.Area, 0.5)
			}
		},
	}

	planimeter29 = planimeterTest{
		testNum:     29,
		description: "Check fix to transitdirect vs transit zero handling inconsistency",
		logic: func(t *testing.T) {
			polygon := NewPolygonArea(WGS84, false)
			polygon.AddPoint(0, 0)
			polygon.AddEdge(90, 1000)
			polygon.AddEdge(0, 1000)
			polygon.AddEdge(-90, 1000)
			a := polygon.Compute(false, true)
			// The area should be 1e6. Prior to the fix it was 1e6 - A/2, where A = ellipsoid area.
			assert.InDelta(t, 1000000.0, a.Area, 0.01)
		},
	}
)
