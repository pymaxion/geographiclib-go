package geodesic

import (
	"fmt"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestConstants(t *testing.T) {
	t.Run("tiny should have expected properties", func(t *testing.T) {
		assert.Greater(t, tiny*epsilon, 0.)
		assert.Equal(t, tiny+epsilon, epsilon)
	})
}

func TestIsFinite(t *testing.T) {
	testCases := []struct {
		val            float64
		expectedResult bool
	}{
		{val: 0, expectedResult: true},
		{val: -0, expectedResult: true},
		{val: 1, expectedResult: true},
		{val: -1, expectedResult: true},
		{val: math.Pi, expectedResult: true},
		{val: math.MaxFloat64 - 1, expectedResult: true},
		{val: math.MaxFloat64, expectedResult: true},
		{val: math.MaxFloat64 + 1, expectedResult: true},
		{val: -math.MaxFloat64 + 1, expectedResult: true},
		{val: -math.MaxFloat64, expectedResult: true},
		{val: -math.MaxFloat64 - 1, expectedResult: true},
		{val: math.Inf(1), expectedResult: false},
		{val: math.Inf(-1), expectedResult: false},
		{val: math.NaN(), expectedResult: false},
	}

	for _, tt := range testCases {
		t.Run(fmt.Sprintf("%.3f", tt.val), func(t *testing.T) {
			assert.Equal(t, tt.expectedResult, isfinite(tt.val))
		})
	}
}
