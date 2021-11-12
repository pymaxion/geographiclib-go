package geodesic

import (
	"testing"
)

func TestPosition(t *testing.T) {
	for _, testCase := range []geodSolve{
		geodSolve17,
		geodSolve80,
	} {
		t.Run(testCase.String(), testCase.logic)
	}
}
