package geodesic

import (
	"fmt"
	"math"
	"testing"

	"github.com/stretchr/testify/assert"
)

var minusZero = math.Copysign(0, -1)

func TestConstants(t *testing.T) {
	t.Run("tiny should have expected properties", func(t *testing.T) {
		assert.Greater(t, tiny*epsilon, 0.)
		assert.Equal(t, tiny+epsilon, epsilon)
	})
}

func TestAngRound(t *testing.T) {
	testCases := []struct {
		val            float64
		expectedResult float64
	}{
		{val: -epsilon / 32, expectedResult: -epsilon / 32},
		{val: -epsilon / 64, expectedResult: minusZero},
		{val: minusZero, expectedResult: minusZero},
		{val: 0.0, expectedResult: 0},
		{val: epsilon / 64, expectedResult: 0},
		{val: epsilon / 32, expectedResult: +epsilon / 32},
		{val: (1 - 2*epsilon) / 64, expectedResult: (1 - 2*epsilon) / 64},
		{val: (1 - epsilon) / 64, expectedResult: 1.0 / 64},
		{val: (1 - epsilon/2) / 64, expectedResult: 1.0 / 64},
		{val: (1 - epsilon/4) / 64, expectedResult: 1.0 / 64},
		{val: 1.0 / 64, expectedResult: 1.0 / 64},
		{val: (1 + epsilon/2) / 64, expectedResult: 1.0 / 64},
		{val: (1 + epsilon) / 64, expectedResult: 1.0 / 64},
		{val: (1 + 2*epsilon) / 64, expectedResult: (1 + 2*epsilon) / 64},
		{val: (1 - epsilon) / 32, expectedResult: (1 - epsilon) / 32},
		{val: (1 - epsilon/2) / 32, expectedResult: 1.0 / 32},
		{val: (1 - epsilon/4) / 32, expectedResult: 1.0 / 32},
		{val: 1.0 / 32, expectedResult: 1.0 / 32},
		{val: (1 + epsilon/2) / 32, expectedResult: 1.0 / 32},
		{val: (1 + epsilon) / 32, expectedResult: (1 + epsilon) / 32},
		{val: (1 - epsilon) / 16, expectedResult: (1 - epsilon) / 16},
		{val: (1 - epsilon/2) / 16, expectedResult: (1 - epsilon/2) / 16},
		{val: (1 - epsilon/4) / 16, expectedResult: 1.0 / 16},
		{val: 1.0 / 16, expectedResult: 1.0 / 16},
		{val: (1 + epsilon/4) / 16, expectedResult: 1.0 / 16},
		{val: (1 + epsilon/2) / 16, expectedResult: 1.0 / 16},
		{val: (1 + epsilon) / 16, expectedResult: (1 + epsilon) / 16},
		{val: (1 - epsilon) / 8, expectedResult: (1 - epsilon) / 8},
		{val: (1 - epsilon/2) / 8, expectedResult: (1 - epsilon/2) / 8},
		{val: (1 - epsilon/4) / 8, expectedResult: 1.0 / 8},
		{val: (1 + epsilon/2) / 8, expectedResult: 1.0 / 8},
		{val: (1 + epsilon) / 8, expectedResult: (1 + epsilon) / 8},
		{val: 1 - epsilon, expectedResult: 1 - epsilon},
		{val: 1 - epsilon/2, expectedResult: 1 - epsilon/2},
		{val: 1 - epsilon/4, expectedResult: 1},
		{val: 1.0, expectedResult: 1},
		{val: 1 + epsilon/4, expectedResult: 1},
		{val: 1 + epsilon/2, expectedResult: 1},
		{val: 1 + epsilon, expectedResult: 1 + epsilon},
		{val: 90.0 - 64*epsilon, expectedResult: 90 - 64*epsilon},
		{val: 90.0 - 32*epsilon, expectedResult: 90},
		{val: 90.0, expectedResult: 90},
	}

	for i, tt := range testCases {
		t.Run(fmt.Sprintf("%d", i), func(t *testing.T) {
			actual := angRound(tt.val)
			assert.True(t, equiv(tt.expectedResult, actual))
		})
	}
}

func TestSincosd(t *testing.T) {
	t.Run("edge cases", func(t *testing.T) {
		type result struct{ sinx, cosx float64 }
		testCases := []struct {
			val            float64
			expectedResult result
		}{
			{val: math.Inf(-1), expectedResult: result{math.NaN(), math.NaN()}},
			{val: -810, expectedResult: result{-1, 0}},
			{val: -720, expectedResult: result{minusZero, 1}},
			{val: -630, expectedResult: result{1, 0}},
			{val: -540, expectedResult: result{minusZero, -1}},
			{val: -450, expectedResult: result{-1, 0}},
			{val: -360, expectedResult: result{minusZero, 1}},
			{val: -270, expectedResult: result{1, 0}},
			{val: -180, expectedResult: result{minusZero, -1}},
			{val: -90, expectedResult: result{-1, 0}},
			{val: minusZero, expectedResult: result{minusZero, 1}},
			{val: 0, expectedResult: result{0, 1}},
			{val: 90, expectedResult: result{1, 0}},
			{val: 180, expectedResult: result{0, -1}},
			{val: 270, expectedResult: result{-1, 0}},
			{val: 360, expectedResult: result{0, 1}},
			{val: 450, expectedResult: result{1, 0}},
			{val: 540, expectedResult: result{0, -1}},
			{val: 630, expectedResult: result{-1, 0}},
			{val: 720, expectedResult: result{0, 1}},
			{val: 810, expectedResult: result{1, 0}},
			{val: math.Inf(1), expectedResult: result{math.NaN(), math.NaN()}},
			{val: math.NaN(), expectedResult: result{math.NaN(), math.NaN()}},
		}

		for _, tt := range testCases {
			t.Run(fmt.Sprintf("%.3f", tt.val), func(t *testing.T) {
				sinx, cosx := sincosd(tt.val)
				assert.True(t, equiv(tt.expectedResult.sinx, sinx))
				assert.True(t, equiv(tt.expectedResult.cosx, cosx))
			})
		}
	})

	t.Run("accuracy", func(t *testing.T) {
		s1, c1 := sincosd(9)
		s2, c2 := sincosd(81)
		s3, c3 := sincosd(-123456789)

		assert.True(t, equiv(s1, c2))
		assert.True(t, equiv(s1, s3))
		assert.True(t, equiv(c1, s2))
		assert.True(t, equiv(c1, -c3))
	})
}

func TestAtan2d(t *testing.T) {
	t.Run("edge cases", func(t *testing.T) {
		testCases := []struct {
			y              float64
			x              float64
			expectedResult float64
		}{
			{0, minusZero, +180},
			{minusZero, minusZero, -180},
			{0, 0, 0},
			{minusZero, 0, minusZero},
			{0, -1, +180},
			{minusZero, -1, -180},
			{0, 1, 0},
			{minusZero, 1, minusZero},
			{-1, 0, -90},
			{-1, minusZero, -90},
			{1, 0, +90},
			{1, minusZero, +90},
			{1, math.Inf(-1), +180},
			{-1, math.Inf(-1), -180},
			{1, math.Inf(1), 0},
			{-1, math.Inf(1), minusZero},
			{math.Inf(1), 1, +90},
			{math.Inf(1), -1, +90},
			{math.Inf(-1), 1, -90},
			{math.Inf(-1), -1, -90},
			{math.Inf(1), math.Inf(-1), +135},
			{math.Inf(-1), math.Inf(-1), -135},
			{math.Inf(1), math.Inf(1), +45},
			{math.Inf(-1), math.Inf(1), -45},
			{math.NaN(), 1, math.NaN()},
			{1, math.NaN(), math.NaN()},
		}

		for _, tt := range testCases {
			t.Run(fmt.Sprintf("y: %.3f, x: %.3f", tt.y, tt.x), func(t *testing.T) {
				assert.True(t, equiv(tt.expectedResult, atan2d(tt.y, tt.x)))
			})
		}
	})

	t.Run("accuracy", func(t *testing.T) {
		s := 7e-16
		assert.Equal(t, atan2d(s, -1), 180-atan2d(s, 1.0))
	})
}

func TestSum(t *testing.T) {
	testCases := []struct {
		u              float64
		v              float64
		expectedResult float64
	}{
		{9, -9, 0},
		{-9, 9, 0},
		{minusZero, 0, 0},
		{0, minusZero, 0},
		{minusZero, minusZero, minusZero},
		{0, 0, 0},
	}

	for i, tt := range testCases {
		t.Run(fmt.Sprintf("edge case %d", i), func(t *testing.T) {
			result, _ := sum(tt.u, tt.v)
			assert.True(t, equiv(tt.expectedResult, result))
		})
	}
}

func TestAngNormalize(t *testing.T) {
	testCases := []struct {
		angle          float64
		expectedResult float64
	}{
		{-900, -180},
		{-720, minusZero},
		{-540, -180},
		{-360, minusZero},
		{-180, -180},
		{minusZero, minusZero},
		{0, 0},
		{180, 180},
		{360, 0},
		{540, 180},
		{720, 0},
		{900, 180},
	}

	for _, tt := range testCases {
		t.Run(fmt.Sprintf("%.3f", tt.angle), func(t *testing.T) {
			assert.True(t, equiv(tt.expectedResult, angNormalize(tt.angle)))
		})
	}
}

func TestAngDiff(t *testing.T) {
	t.Run("edge cases", func(t *testing.T) {
		testCases := []struct {
			x              float64
			y              float64
			expectedResult float64
		}{
			{0, 0, 0},
			{0, minusZero, minusZero},
			{minusZero, 0, 0},
			{minusZero, minusZero, 0},
			{5, 365, 0},
			{365, 5, minusZero},
			{5, 185, 180},
			{185, 5, -180},
			{epsilon, 180, 180},
			{-epsilon, 180, -180},
			{epsilon, -180, 180},
			{-epsilon, -180, -180},
		}

		for _, tt := range testCases {
			t.Run(fmt.Sprintf("x: %.3f, y: %.3f", tt.x, tt.y), func(t *testing.T) {
				d, _ := angDiff(tt.x, tt.y)
				assert.True(t, equiv(tt.expectedResult, d))
			})
		}
	})

	t.Run("accuracy", func(t *testing.T) {
		x := 138 + 128*epsilon
		y := -164.
		d, _ := angDiff(x, y)
		assert.Equal(t, 58-128*epsilon, d)
	})
}

func TestIsFinite(t *testing.T) {
	testCases := []struct {
		val            float64
		expectedResult bool
	}{
		{val: 0, expectedResult: true},
		{val: minusZero, expectedResult: true},
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

// equiv is a helper function to test for equivalence
func equiv(x, y float64) bool {
	return (math.IsNaN(x) && math.IsNaN(y)) || (x == y && math.Copysign(1, x) == math.Copysign(1, y))
}
