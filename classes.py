import re


class Integer:
    """
    Класс для работы с целыми числами, представленными в виде строки.
    Хранит число как массив цифр и его знак.
    """

    def __init__(self, number: str):
        """
        Инициализация объекта целого числа.

        :param number: Строка, представляющая целое число.
        """
        if not isinstance(number, str) or not number.lstrip('-').isdigit():
            raise ValueError("Некорректный формат: требуется строка, представляющая целое число.")

        self.is_negative = number.startswith('-')
        self.digits = [int(d) for d in number.lstrip('-')]

    def __str__(self):
        """
        Возвращает строковое представление числа.
        """
        sign = '-' if self.is_negative else ''
        return sign + ''.join(map(str, self.digits))

    def __int__(self):
        """
        Возвращает число в виде целого типа Python.
        """
        return int(str(self))

    def get_digits(self):
        """Возвращает массив цифр числа."""
        return self.digits

    def get_sign(self):
        """Возвращает знак числа (True для отрицательных)."""
        return self.is_negative

    def __len__(self):
        """Возвращает количество цифр в числе."""
        return len(self.digits)

    def __eq__(self, other):
        if isinstance(other, Integer):
            return str(self) == str(other)
        return False


class Natural:
    """
    Класс для работы с натуральными числами и нулем.
    """

    def __init__(self, number: str):
        """
        Инициализация натурального числа.

        :param number: Строка, представляющая натуральное число.
        """
        if not isinstance(number, str) or not number.isdigit() or int(number) < 0:
            raise ValueError("Некорректный формат: требуется строка, представляющая натуральное число или ноль.")

        self.digits = [int(d) for d in str(int(number))]

    def __str__(self):
        """
        Возвращает строковое представление числа.
        """
        return ''.join(map(str, self.digits))

    def __int__(self):
        """
        Возвращает число в виде целого типа Python.
        """
        return int(str(self))

    def get_digits(self):
        """Возвращает массив цифр числа."""
        return self.digits

    def __len__(self):
        """Возвращает количество цифр в числе."""
        return len(self.digits)

    def __eq__(self, other):
        if isinstance(other, Natural):
            return str(self) == str(other)
        return False


class Rational:
    """
    Класс для работы с рациональными числами.
    Хранит числитель как Integer и знаменатель как Natural.
    """

    def __init__(self, numerator: Integer, denominator: Natural = Natural("1")):
        """
        Инициализация рационального числа.

        :param numerator: Числитель (Integer).
        :param denominator: Знаменатель (Natural), по умолчанию равен 1.
        """
        if int(denominator) == 0:
            raise ValueError("Знаменатель не может быть нулём.")

        self.numerator = numerator
        self.denominator = denominator

    def __str__(self):
        """
        Возвращает строковое представление рационального числа.
        """
        return f"{self.numerator}/{self.denominator}" if int(self.denominator) != 1 else str(self.numerator)

    def __eq__(self, other):
        if isinstance(other, Rational):
            return self.numerator == other.numerator and self.denominator == other.denominator
        return False


class Polynomial:
    """
    Класс для работы с многочленами на основе словаря.
    """

    def __init__(self):
        """Инициализация пустого многочлена"""
        self.terms = {}  # {degree (int): coefficient (Rational)}
        self._degree = None
        self._leading_coeff = None

    def add_term(self, degree: Natural, coefficient: Rational):
        """Добавление/обновление члена многочлена"""
        deg = int(degree)
        if int(coefficient.numerator) == 0:
            self.terms.pop(deg, None)
        else:
            self.terms[deg] = coefficient

        # Инвалидация кэша
        self._degree = None
        self._leading_coeff = None

    def __str__(self):
        """Строковое представление многочлена"""
        if not self.terms:
            return "0"

        result = []
        for degree in sorted(self.terms.keys(), reverse=True):
            coeff = self.terms[degree]

            if int(coeff.numerator) == 0:
                continue

            sign = " + " if int(coeff.numerator) > 0 and result else ""
            if int(coeff.numerator) < 0:
                sign = " - " if result else "-"

            abs_num = abs(int(coeff.numerator))
            if int(coeff.denominator) != 1:
                coeff_str = f"{abs_num}/{coeff.denominator}"
            else:
                coeff_str = str(abs_num)

            if degree > 1:
                term = f"{coeff_str}x^{degree}"
            elif degree == 1:
                term = f"{coeff_str}x"
            else:
                term = coeff_str

            result.append(sign + term)

        return "".join(result)

    def build_from_string(self, input_str: str):
        """Построение многочлена из строки"""
        input_str = input_str.replace(' ', '').replace('*', '')
        input_str = input_str.replace('(', '').replace(')', '')

        terms = re.split(r'(?=[+-])', input_str)
        for term in terms:
            if not term:
                continue
            term = term.lstrip('+')

            if 'x' in term:
                if '^' in term:
                    coeff, _, degree = term.partition('x^')
                else:
                    coeff, _, degree = term.partition('x')
                    degree = "1"
                coeff = coeff or "1"
                if coeff == "-":
                    coeff = "-1"
            else:
                coeff = term
                degree = "0"

            try:
                if '/' in coeff:
                    num, den = coeff.split('/')
                    is_negative = num.startswith('-')
                    num = num.lstrip('-')
                    numerator = -int(num) if is_negative else int(num)
                    denominator = int(den)
                else:
                    numerator = int(coeff)
                    denominator = 1

                rational_coeff = Rational(Integer(str(numerator)), Natural(str(denominator)))
                self.add_term(Natural(str(degree)), rational_coeff)

            except ValueError as e:
                raise ValueError(f"Ошибка при разборе коэффициента '{coeff}': {str(e)}")

    def getCoeff(self, deg: Natural) -> Rational:
        """Получение коэффициента при заданной степени"""
        degree = int(deg)
        return self.terms.get(degree, Rational(Integer("0"), Natural("1")))

    def get_degree(self) -> Natural:
        """Получение степени многочлена"""
        if self._degree is None:
            self._degree = max(self.terms.keys()) if self.terms else 0
        return Natural(str(self._degree))

    def get_leading_coeff(self) -> Rational:
        """Получение старшего коэффициента"""
        if self._leading_coeff is None:
            deg = self.get_degree()
            self._leading_coeff = self.terms.get(int(deg), Rational(Integer("0")))
        return self._leading_coeff

    def __eq__(self, other):
        if not isinstance(other, Polynomial):
            return False
        return str(self) == str(other)


def create_polynomial(input_str: str) -> Polynomial:
    """Утилита для создания многочлена из строки"""
    poly = Polynomial()
    if input_str == "0":
        return poly
    poly.build_from_string(input_str)
    return poly