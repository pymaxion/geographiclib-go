package inverse

import (
	"bufio"
	"fmt"
	"os"

	"geographiclib-go/examples/util"
	"geographiclib-go/geodesic"
)

func main() {
	var lat1, lon1, lat2, lon2 float64
	var err error
	scanner := bufio.NewScanner(os.Stdin)

	for {
		fmt.Println()
		fmt.Println("Solve the inverse geodesic problem. Enter values lat1, lon1, lat2, lon2:")
		if lat1, err = util.PromptFloat(scanner, "lat1"); err != nil {
			fmt.Println(err.Error())
			continue
		}
		if lon1, err = util.PromptFloat(scanner, "lon1"); err != nil {
			fmt.Println(err.Error())
			continue
		}
		if lat2, err = util.PromptFloat(scanner, "lat2"); err != nil {
			fmt.Println(err.Error())
			continue
		}
		if lon2, err = util.PromptFloat(scanner, "lon2"); err != nil {
			fmt.Println(err.Error())
			continue
		}

		r := geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2)
		fmt.Printf("azi1: %.6f, azi2: %.6f, s12: %.6f", r.Azi1, r.Azi2, r.S12)
		fmt.Println()
	}
}
