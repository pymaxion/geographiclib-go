package geodesic

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestConstants(t *testing.T) {
	t.Run("tiny should have expected properties", func(t *testing.T) {
		assert.Greater(t, tiny*epsilon, 0.)
		assert.Equal(t, tiny+epsilon, epsilon)
	})
}
