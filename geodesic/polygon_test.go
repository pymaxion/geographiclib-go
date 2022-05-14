package geodesic

import "testing"

func TestPolygonArea(t *testing.T) {
	for _, testCase := range []planimeterTest{
		planimeter0,
		planimeter5,
		planimeter6,
		planimeter12,
		planimeter12r,
		planimeter13,
		planimeter15,
		planimeter19,
		planimeter21,
		planimeter29,
	} {
		t.Run(testCase.String(), testCase.logic)
	}
}
