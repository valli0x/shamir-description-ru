package main

import (
	"crypto/aes"
	"crypto/rand"
	"fmt"
	"io"
	"shamir-description-ru/shamir"
)

func main() {
	splitKey, err := GenerateKey(rand.Reader)
	if err != nil {
		fmt.Println(err)
	}
	fmt.Println("splitKey:\n", splitKey)

	keys, err := shamir.Split(splitKey, 3, 2)
	if err != nil {
		fmt.Println(err)
	}

	fmt.Println("keys:")
	for _, key := range keys {
		fmt.Println(key)
	}
}

func GenerateKey(reader io.Reader) ([]byte, error) {
	// Generate a 256bit key
	buf := make([]byte, 2*aes.BlockSize)
	_, err := reader.Read(buf)

	return buf, err
}
