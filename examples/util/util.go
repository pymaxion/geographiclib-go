package util

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"strconv"
)

func PromptFloat(scanner *bufio.Scanner, name string) (float64, error) {
	fmt.Print(name + ": ")
	return nextFloat(scanner)
}

func nextFloat(scanner *bufio.Scanner) (float64, error) {
	if !scanner.Scan() {
		fmt.Println("Exiting...")
		if scanner.Err() != nil {
			log.Fatal(scanner.Err().Error())
		}
		os.Exit(0)
	}
	return strconv.ParseFloat(scanner.Text(), 64)
}
