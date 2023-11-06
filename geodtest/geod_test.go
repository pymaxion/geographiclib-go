package main

import (
	"bufio"
	"errors"
	"fmt"
	"os"
	"strconv"
	"strings"
	"testing"

	"github.com/pymaxion/geographiclib-go/geodesic"
	"github.com/stretchr/testify/assert"
)

const geodTestDatFilepathEnvVar = "GEODTEST_DAT_PATH"

func Test_GeodTest(t *testing.T) {
	geodTestDatFilepath := os.Getenv(geodTestDatFilepathEnvVar)
	if geodTestDatFilepath == "" {
		t.Skipf("\n"+
			"Skipped full GeodTest test suite because %s was not set. To run this test suite, perform the following steps:\n"+
			"  1. Download GeodTest.dat.gz from https://sourceforge.net/projects/geographiclib/files/testdata/GeodTest.dat.gz\n"+
			"  2. gunzip GeodTest.dat.gz to produce GeodTest.dat\n"+
			"  3. Set the %s environment variable to the filepath of GeodTest.dat on your machine\n"+
			"Once %s is set to a valid filepath, please rerun the tests.\n",
			geodTestDatFilepathEnvVar, geodTestDatFilepathEnvVar, geodTestDatFilepathEnvVar)
	}

	f, err := os.Open(geodTestDatFilepath)
	if err != nil {
		t.Fatal(err)
	}
	fs, err := f.Stat()
	if err != nil {
		t.Fatal(err)
	}
	if fs.IsDir() {
		t.Fatalf("cannot read %s as it is not a file", f.Name())
	}

	scanner := bufio.NewScanner(f)
	scanner.Split(bufio.ScanLines)

	for i := 0; scanner.Scan(); i++ {
		line := scanner.Text()
		ps, err := parse(line)
		if err != nil {
			t.Fatal(err)
		}

		t.Run(fmt.Sprintf("GeodTest.dat #%d", i), func(t *testing.T) {
			t.Run("direct (from point 1)", func(t *testing.T) {
				t.Parallel()
				r := geodesic.WGS84.Direct(ps.lat1, ps.lon1, ps.azi1, ps.s12)

				assert.InDelta(t, ps.lat2, r.Lat2, 5e-6)
				assert.InDelta(t, ps.lon2, r.Lon2, 5e-6)
				assert.InDelta(t, ps.azi2, r.Azi2, 5e-6)
				assert.InDelta(t, ps.a12, r.A12, 5e-6)
			})

			t.Run("direct (from point 2)", func(t *testing.T) {
				t.Parallel()
				r := geodesic.WGS84.Direct(ps.lat2, ps.lon2, ps.azi2, -ps.s12)

				assert.InDelta(t, ps.lat1, r.Lat2, 5e-6)
				assert.InDelta(t, ps.lon1, r.Lon2, 5e-6)
				assert.InDelta(t, ps.azi1, r.Azi2, 5e-6)
				assert.InDelta(t, -ps.a12, r.A12, 5e-6)
			})

			t.Run("inverse", func(t *testing.T) {
				t.Parallel()
				r := geodesic.WGS84.Inverse(ps.lat1, ps.lon1, ps.lat2, ps.lon2)

				if ps.azi1 > 89.5 {
					assert.InDelta(t, ps.azi1, r.Azi1, 0.02)
				} else {
					assert.InDelta(t, ps.azi1, r.Azi1, 1e-5)
				}

				if ps.azi2 > 89.5 {
					assert.InDelta(t, ps.azi2, r.Azi2, 0.02)
				} else {
					assert.InDelta(t, ps.azi2, r.Azi2, 1e-5)
				}

				assert.InDelta(t, ps.s12, r.S12, 5e-6)
				assert.InDelta(t, ps.a12, r.A12, 5e-6)
			})
		})
	}

	err = scanner.Err()
	if err != nil {
		t.Fatal(err)
	}
}

type params struct {
	lat1       float64
	lon1       float64
	azi1       float64
	lat2       float64
	lon2       float64
	azi2       float64
	s12        float64
	a12        float64
	m12Reduced float64
	s12Area    float64
}

func parse(line string) (params, error) {
	pieces := strings.Split(line, " ")
	if len(pieces) != 10 {
		return params{}, errors.New(fmt.Sprintf("unexpectedly formatted geodtest line: '%s'", line))
	}

	vals := make([]float64, len(pieces))
	var err error
	for i := 0; i < len(pieces); i++ {
		vals[i], err = strconv.ParseFloat(pieces[i], 64)
		if err != nil {
			return params{}, err
		}
	}

	return params{
		lat1:       vals[0],
		lon1:       vals[1],
		azi1:       vals[2],
		lat2:       vals[3],
		lon2:       vals[4],
		azi2:       vals[5],
		s12:        vals[6],
		a12:        vals[7],
		m12Reduced: vals[8],
		s12Area:    vals[9],
	}, nil
}
