package main

import (
	"bufio"
	"fmt"
	"os"

	"github.com/pymaxion/geographiclib-go/examples/util"
	"github.com/pymaxion/geographiclib-go/geodesic"
)

func main() {
	var lat1, lon1, azi1, s12 float64
	var err error
	scanner := bufio.NewScanner(os.Stdin)

	for {
		fmt.Println()
		fmt.Println("Solve the direct geodesic problem. Enter values lat1, lon1, azi1, s12:")
		if lat1, err = util.PromptFloat(scanner, "lat1"); err != nil {
			fmt.Println(err.Error())
			continue
		}
		if lon1, err = util.PromptFloat(scanner, "lon1"); err != nil {
			fmt.Println(err.Error())
			continue
		}
		if azi1, err = util.PromptFloat(scanner, "azi1"); err != nil {
			fmt.Println(err.Error())
			continue
		}
		if s12, err = util.PromptFloat(scanner, "s12"); err != nil {
			fmt.Println(err.Error())
			continue
		}

		r := geodesic.WGS84.Direct(lat1, lon1, azi1, s12)
		fmt.Printf("lat2: %.6f, lon2: %.6f, azi2: %.6f", r.Lat2, r.Lon2, r.Azi2)
		fmt.Println()
	}
}
