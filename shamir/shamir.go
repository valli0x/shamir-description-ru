package shamir

import (
	"crypto/rand"
	"crypto/subtle"
	"fmt"
	mathrand "math/rand"
	"time"
)

const (
	// ShareOverhead is the byte size overhead of each share
	// when using Split on a secret. This is caused by appending
	// a one byte tag to the share.
	ShareOverhead = 1
)

// polynomial represents a polynomial of arbitrary degree
type polynomial struct {
	coefficients []uint8
}

/* 
	Здесь все просто, степень полинома нам известна - это degree
	все полиномы выглядят одинаково
	-- формула или картинка полинома --
	те нам нужно заполнить случайными коэффициентами у каждого x

	создаем слайс с размером degree+1, тк нужно сохранить наш байт(он же свободный член)
	и заполняем его случайнымы байтами, кроме нашего байта или свободного члена в полиноме
*/

// makePolynomial constructs a random polynomial of the given
// degree but with the provided intercept value.
func makePolynomial(intercept, degree uint8) (polynomial, error) {
	// Create a wrapper
	p := polynomial{
		coefficients: make([]byte, degree+1),
	}

	// Ensure the intercept is set
	p.coefficients[0] = intercept

	// Assign random co-efficients to the polynomial
	if _, err := rand.Read(p.coefficients[1:]); err != nil {
		return p, err
	}

	return p, nil
}

/*
	Нам передают x и нам нужно для этого полинома посчитать значение y
	В лоб при больших степенях считать сложно. Есть инструмент - схема Горнера
	Подробнее о ней - https://youtu.be/H9teIpb62iI?t=180

	кратко: есть наш полином 4 степени - 2x^4 + 7x^3 - 2x^2 - 13x + 6
	нам x умножать на коэффициент с начала уравнения и прибавлять следующий
	т.е. x=1 
	1) 1*2+7=9
	2) 1*9+(-2)=7
	3) 1*7+(-13)=-6
	4) 1*(-6)+6=0
	наш y=0(он может быть и не равен 0, с видео взял пример, там корни уравнения искали)

	соответственно при x=0 у нас получится наш свободный член, поэтому не считая сразу его вернем
	мы берем последний элемент(в нашем слайсе коэффициенты расположены с большого индекса к меньшему
	соответствуя коэффициентам у многочлена от большей степени к меньшей) и умнажая на x, прибавляя следующий
	коэффициент
	так проходимся по всем коэффициентам, умножая x на результат и прибавляя следующий, пока не получим наш y

	ВАЖНАЯ ЧАСТЬ:
	функции add(a, b uint8) mult(a, b uint8) div(a, b uint8) используются для сложения вычитания умножения и деления
	add() - для сложения и вычитания вместе
	Мы не должны выходить за пределы байта, те 256 значений, поэтому все операции происходят в так называемом конечном поле
	Это довольно долгая история, поэтому я решил просто условится на том, что у нас есть эти функции
	Сами hashcorp предоставляют эту ссылку - http://www.samiam.org/galois.html если интересно устройство этих операций
*/

// evaluate returns the value of the polynomial for the given x
func (p *polynomial) evaluate(x uint8) uint8 {
	// Special case the origin
	if x == 0 {
		return p.coefficients[0]
	}

	// Compute the polynomial value using Horner's method.
	degree := len(p.coefficients) - 1
	out := p.coefficients[degree]
	for i := degree - 1; i >= 0; i-- {
		coeff := p.coefficients[i]
		out = add(mult(out, x), coeff)
	}
	return out
}

/*
	первый и второй аргументы - это набор x и соответствующий им y, 3 - это x для которого
	нужно найти y в полиноме

	Подробно об полиноме Лагранжа - https://www.youtube.com/watch?v=nj2RBZ6xFeY
	Не пугайтесь базисного полинома, если уделить минут 15 для разбора, все окажется легче
	Наша формула:
	--- скрин из видео ---

	Нужно последовательно умножить все y на так называемый базисный полином и сложить результаты
	Первый цикл - занимается умножением y на базисный полином и сложением - т.е. отражение формулы внизу
	group - умножение y на базисный полином
	result - результат сложения

	Второй(вложенный) - высчитывет базисный полином - это формулы сверху.
	Еще раз напомню:
	add(a, b uint8) - это сложение и вычитание в конечном поле, они одинаковы
	div(a, b uint8) - деление, mult(a, b uint8) - умножение

	Базисный полином для y это несколько перемножаемых дробей:
	num - результат вычисления числителя для n дроби
	(во всех дробях x минус последовательно перебираемые значения точек x, пропуская тот x который соответствует y, условие i==j)
	в нашем случае x известен, мы хотим получить результат для x=0, поэтому 0 - перебираемые значения...

	denom - результат вычисления знаменателя для n дроби
	(перебираемый x минус последовательно итерируемые значения x, пропуская тот x который соответствует y, условие i==j)

	примерно так выглядит для значения x под индексом 1 в слайсе:
	(x-x[0])/(x[1]-x[0]) * (x-x[2])/(x[1]-x[2]) * (x-x[3])/(x[1]-x[3])

	пропустили (x-x[1])/(x[1]-x[1]) по правилу, да и ноль в знаменателе получается к тому же
	поскольку нужен результат полинома при x=0, чтобы получить свободный член, то базисный полином получается
	(0-x[0])/(x[1]-x[0]) * (0-x[2])/(x[1]-x[2]) * (0-x[3])/(x[1]-x[3])

	term - результат деления num на denom, т.е. результат одной дроби(к примеру (x-x[0])/(x[1]-x[0]))
	basis - все результаты дробей(term) базисного полинома перемножаются
*/

// interpolatePolynomial takes N sample points and returns
// the value at a given x using a lagrange interpolation.
func interpolatePolynomial(x_samples, y_samples []uint8, x uint8) uint8 {
	limit := len(x_samples)
	var result, basis uint8
	for i := 0; i < limit; i++ {
		basis = 1
		for j := 0; j < limit; j++ {
			if i == j {
				continue
			}
			num := add(x, x_samples[j])
			denom := add(x_samples[i], x_samples[j])
			term := div(num, denom)
			basis = mult(basis, term)
		}
		group := mult(y_samples[i], basis)
		result = add(result, group)
	}
	return result
}

// div divides two numbers in GF(2^8)
func div(a, b uint8) uint8 {
	if b == 0 {
		// leaks some timing information but we don't care anyways as this
		// should never happen, hence the panic
		panic("divide by zero")
	}

	log_a := logTable[a]
	log_b := logTable[b]
	diff := ((int(log_a) - int(log_b)) + 255) % 255

	ret := int(expTable[diff])

	// Ensure we return zero if a is zero but aren't subject to timing attacks
	ret = subtle.ConstantTimeSelect(subtle.ConstantTimeByteEq(a, 0), 0, ret)
	return uint8(ret)
}

// mult multiplies two numbers in GF(2^8)
func mult(a, b uint8) (out uint8) {
	log_a := logTable[a]
	log_b := logTable[b]
	sum := (int(log_a) + int(log_b)) % 255

	ret := int(expTable[sum])

	// Ensure we return zero if either a or b are zero but aren't subject to
	// timing attacks
	ret = subtle.ConstantTimeSelect(subtle.ConstantTimeByteEq(a, 0), 0, ret)
	ret = subtle.ConstantTimeSelect(subtle.ConstantTimeByteEq(b, 0), 0, ret)

	return uint8(ret)
}

// add combines two numbers in GF(2^8)
// This can also be used for subtraction since it is symmetric.
func add(a, b uint8) uint8 {
	return a ^ b
}

// secret - это наш ключ представленный в виде байтов, parts - общее кол-во частей
// threshold - пороговое кол-во частей, возвращаемые значения - это наши ключи и ошибка

// Split takes an arbitrarily long secret and generates a `parts`
// number of shares, `threshold` of which are required to reconstruct
// the secret. The parts and threshold must be at least 2, and less
// than 256. The returned shares are each one byte longer than the secret
// as they attach a tag used to reconstruct the secret.
func Split(secret []byte, parts, threshold int) ([][]byte, error) {
	
	// Split мы делим на две части: Валидация, подготовка и сама разбивка ключа на несколько частей

	/*
		Валидация:
		1.  Общее кол-во частей не должно быть меньше порогового
		2.  Кол-во частей не должно быть больше 255
			--
			У нас есть только один байт, т.е. 256 комбинаций(от 0 до 255)
			то последующие ключи будут повторяться, т.к. при одинаковых x на конце ключей,
			все значения y под каждым многочленом одинаковы
			--
		3.  Порог не должен быть меньше 2
		4.  Порог не должен быть больше 255
		5.  Ключ не должен быть пустым(они используют имя переменной для ключа - secret,
			у них все секретное - секрет)
	*/

	// Sanity check the input
	if parts < threshold {
		return nil, fmt.Errorf("parts cannot be less than threshold")
	}
	if parts > 255 {
		return nil, fmt.Errorf("parts cannot exceed 255")
	}
	if threshold < 2 {
		return nil, fmt.Errorf("threshold must be at least 2")
	}
	if threshold > 255 {
		return nil, fmt.Errorf("threshold cannot exceed 255")
	}
	if len(secret) == 0 {
		return nil, fmt.Errorf("cannot split an empty secret")
	}

	/*
		Создаем слайс(расширяющийся массив в golang) из случайных x 
		до этого просто брали от 1 до 5, но потом заменили на случайные значения
		https://github.com/hashicorp/vault/commit/b4602fc2443dc3de451cd9aeac7c4d79390adbd8

		И слайсы для каждого ключа, размером с ключ + 1 байт для x на конце, и ставим x из случайного набора
		в конец слайса
	*/

	// Generate random list of x coordinates
	mathrand.Seed(time.Now().UnixNano())
	xCoordinates := mathrand.Perm(255)

	// Allocate the output array, initialize the final byte
	// of the output with the offset. The representation of each
	// output is {y1, y2, .., yN, x}.
	out := make([][]byte, parts)
	for idx := range out {
		out[idx] = make([]byte, len(secret)+1)
		out[idx][len(secret)] = uint8(xCoordinates[idx]) + 1
	}


	/*
		Вторая часть это проход по каждому байту ключа, создание полинома со степенью порог-1 для этого байта,
		где наш байт это свободный член.

		Далее берем наши x с концов ключей(из того же массива случайных x), и вычисляем для каждого ключа его
		y. Один ключ это набор из одного x на конце и множества y для каждого построенного полинома для каждого байта
		искомого ключа.

		1. Строим полином для байта - передали байт и степень полинома(порог-1) в функцию makePolynomial
		2. Вычисляем для каждого ключа свой байт при заданном для этого ключа x(он один для ключа) evaluate функция
	*/

	// -- картинка с разбиением --

	// Construct a random polynomial for each byte of the secret.
	// Because we are using a field of size 256, we can only represent
	// a single byte as the intercept of the polynomial, so we must
	// use a new polynomial for each byte.
	for idx, val := range secret {
		p, err := makePolynomial(val, uint8(threshold-1))
		if err != nil {
			return nil, fmt.Errorf("failed to generate polynomial: %w", err)
		}

		// Generate a `parts` number of (x,y) pairs
		// We cheat by encoding the x value once as the final index,
		// so that it only needs to be stored once.
		for i := 0; i < parts; i++ {
			x := uint8(xCoordinates[i]) + 1
			y := p.evaluate(x)
			out[i][idx] = y
		}
	}

	// Return the encoded secrets
	return out, nil
}

/*
	parts - это ключи для объединения
	возвращаемые значения - это исходный ключ и ошибка
*/

// Combine is used to reverse a Split and reconstruct a secret
// once a `threshold` number of parts are available.
func Combine(parts [][]byte) ([]byte, error) {


	/*
		Combine - также мы можем разделить на 2 части
		1. Валидаия и подготовка дынных
		2. Нахождения исходного ключа
	*/

	/*
		Валидация:
		1. Кол-во частей не должно быть меньше двух
		2. Первый ключ в списке не должен быть меньше 2 байт(y,x - минимально 2 точки нужны) 
		и все ключи должны иметь равную длину(сравнивают все с первым)
	*/

	// Verify enough parts provided
	if len(parts) < 2 {
		return nil, fmt.Errorf("less than two parts cannot be used to reconstruct the secret")
	}

	// Verify the parts are all the same length
	firstPartLen := len(parts[0])
	if firstPartLen < 2 {
		return nil, fmt.Errorf("parts must be at least two bytes")
	}
	for i := 1; i < len(parts); i++ {
		if len(parts[i]) != firstPartLen {
			return nil, fmt.Errorf("all parts must be the same length")
		}
	}

	/* 
		secret - наш исходный ключ, поскольку у нас x на конце каждого ключа, иходный будет на один байт меньше
		далее строится список из точек(для x один для y второй)
		ключи(например):
		[23,7,2,45...5] 
		[32,2,3,21...2]
		[18,9,32,1...15]
		берем с концов x - 5, 2, 15 и заполняем для x слайс
		(проверям чтобы x не совпадали иначе ключи тоже будут одинаковы, 
		потому что в полиноме одному x соответсвует один результат он же y)
	*/

	// Create a buffer to store the reconstructed secret
	secret := make([]byte, firstPartLen-1)

	// Buffer to store the samples
	x_samples := make([]uint8, len(parts))
	y_samples := make([]uint8, len(parts))

	// Set the x value for each sample and ensure no x_sample values are the same,
	// otherwise div() can be unhappy
	checkMap := map[byte]bool{}
	for i, part := range parts {
		samp := part[firstPartLen-1]
		if exists := checkMap[samp]; exists {
			return nil, fmt.Errorf("duplicate part detected")
		}
		checkMap[samp] = true
		x_samples[i] = samp
	}

	/*
		Нам нужно восстановить каждый байт исходного ключа, набор из x  у нас есть 
		Мы для первого байта заполняем y
		ключи(например):
		[23,7,2,45...5] 
		[32,2,3,21...2]
		[18,9,32,1...15]

		Для первого байта нам нужен первый столбец, получаются точки - [5,23],[2,32],[15,18]
		Для второго - второй столбец [5,7],[2,2],[15,9] и т.д.
		Далее мы прокидываем в фукнцию interpolatePolynomial наш список из x и y 
		3 - аргумент - это точка для интерполяции, те x для которого нужно найти y в нашем многочлене
		нам нужен свободный коэффициент - наш байт, а он получается при x = 0
	*/

	// Reconstruct each byte
	for idx := range secret {
		// Set the y value for each sample
		for i, part := range parts {
			y_samples[i] = part[idx]
		}

		// Interpolate the polynomial and compute the value at 0
		val := interpolatePolynomial(x_samples, y_samples, 0)

		// Evaluate the 0th value to get the intercept
		secret[idx] = val
	}
	return secret, nil
}
