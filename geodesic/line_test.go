package geodesic

import (
	"testing"
)

func TestLine(t *testing.T) {
	for _, testCase := range []geodSolve{
		geodSolve17,
		geodSolve61,
		geodSolve65,
		geodSolve69,
		geodSolve71,
		geodSolve80,
	} {
		t.Run(testCase.String(), testCase.logic)
	}
}
