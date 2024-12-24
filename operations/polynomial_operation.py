from classes import *
from operations.rational_operations import RationalOperations
from operations.integer_operations import IntegerOperations
from operations.natural_operations import NaturalOperations

class PolynomialOperations:
    @staticmethod
    def ADD_PP_P(poly1: Polynomial, poly2: Polynomial) -> Polynomial:
        """
        Сложение двух многочленов.
        """
        result = Polynomial()

        # Добавляем коэффициенты из poly1
        for degree, coeff in poly1.terms.items():
            result.terms[degree] = coeff

        # Добавляем коэффициенты из poly2
        for degree, coeff in poly2.terms.items():
            if degree in result.terms:
                result.terms[degree] = RationalOperations.ADD_QQ_Q(result.terms[degree], coeff)
            else:
                result.terms[degree] = coeff

        return result

    @staticmethod
    def SUB_PP_P(poly1: Polynomial, poly2: Polynomial) -> Polynomial:
        """
        Вычитание двух многочленов.
        """
        result = Polynomial()

        # Добавляем коэффициенты из poly1
        for degree, coeff in poly1.terms.items():
            result.terms[degree] = coeff

        # Вычитаем коэффициенты из poly2
        for degree, coeff in poly2.terms.items():
            if degree in result.terms:
                result.terms[degree] = RationalOperations.SUB_QQ_Q(result.terms[degree], coeff)
            else:
                result.terms[degree] = RationalOperations.MUL_QQ_Q(coeff, Rational(Integer("-1"), Natural("1")))

        return result

    @staticmethod
    def MUL_PQ_P(poly: Polynomial, q: Rational) -> Polynomial:
        """
        Умножение многочлена на рациональное число.
        """
        result = Polynomial()

        for degree, coeff in poly.terms.items():
            result.terms[degree] = RationalOperations.MUL_QQ_Q(coeff, q)

        return result

    @staticmethod
    def MUL_Pxk_P(poly: Polynomial, k: Natural) -> Polynomial:
        """
        Умножение многочлена на x^k.
        """
        result = Polynomial()

        for degree, coeff in poly.terms.items():
            result.terms[degree + int(k)] = coeff

        return result

    @staticmethod
    def LED_P_Q(poly: Polynomial) -> Rational:
        """
        Возвращает старший коэффициент многочлена.
        """
        max_degree = max(poly.terms.keys())
        return poly.terms[max_degree]

    @staticmethod
    def DEG_P_N(poly: Polynomial) -> int:
        """
        Возвращает степень многочлена.
        """
        return max(poly.terms.keys())

    @staticmethod
    def FAC_P_Q(poly: Polynomial) -> (Rational, Rational):
        """
        Находит НОД числителей и НОК знаменателей коэффициентов многочлена.
        """
        numerators = []
        denominators = []

        # Извлекаем числители и знаменатели всех коэффициентов
        for coeff in poly.terms.values():
            numerators.append(coeff.numerator)
            denominators.append(coeff.denominator)

        # Вычисляем НОД числителей
        gcf_numerators = numerators[0]
        for num in numerators[1:]:
            gcf_numerators = NaturalOperations.GCF_NN_N(gcf_numerators, num)

        # Вычисляем НОК знаменателей
        lcm_denominators = denominators[0]
        for denom in denominators[1:]:
            lcm_denominators = NaturalOperations.LCM_NN_N(lcm_denominators, denom)

        # Возвращаем НОД числителей и НОК знаменателей как объекты Rational
        return Rational(Integer(gcf_numerators), Natural(1)), Rational(Integer(lcm_denominators),Natural(1))

    @staticmethod
    def MUL_PP_P(poly1: Polynomial, poly2: Polynomial) -> Polynomial:
        """
        Умножение двух многочленов.
        """
        result = Polynomial()

        for deg1, coeff1 in poly1.terms.items():
            for deg2, coeff2 in poly2.terms.items():
                new_degree = deg1 + deg2
                new_coeff = RationalOperations.MUL_QQ_Q(coeff1, coeff2)
                if new_degree in result.terms:
                    result.terms[new_degree] = RationalOperations.ADD_QQ_Q(result.terms[new_degree], new_coeff)
                else:
                    result.terms[new_degree] = new_coeff

        return result

    @staticmethod
    def DIV_PP_P(poly1: Polynomial, poly2: Polynomial) -> Polynomial:
        """
        Деление двух многочленов poly1 и poly2. Возвращает только частное.
        """
        # Копирование полинома для остатка
        remainder = Polynomial()
        remainder.terms = poly1.terms.copy()

        # Создаем пустой полином для частного
        quotient = Polynomial()

        # Пока степень остатка больше или равна степени делителя
        while remainder.get_degree() >= poly2.get_degree():
            # Определяем степень и коэффициент для нового члена частного
            degree_diff = remainder.get_degree() - poly2.get_degree()
            coeff = remainder.getCoeff(remainder.get_degree()) / poly2.getCoeff(poly2.get_degree())

            # Новый член частного
            new_term = Polynomial()
            new_term.add_term(degree_diff, coeff)
            quotient.add_term(degree_diff, coeff)

            # Вычитаем из остатка результат умножения poly2 на найденный коэффициент
            temp = PolynomialOperations.MUL_PP_P(poly2, new_term)
            remainder = remainder - temp

        # Возвращаем только частное
        return quotient

    @staticmethod
    def MOD_PP_P(dividend: Polynomial, divisor: Polynomial) -> Polynomial:
        """
        Остаток от деления многочлена на многочлен.
        """
        quotient = PolynomialOperations.DIV_PP_P(dividend, divisor)
        product = PolynomialOperations.MUL_PP_P(quotient, divisor)
        return PolynomialOperations.SUB_PP_P(dividend, product)

    @staticmethod
    def GCF_PP_P(poly1: Polynomial, poly2: Polynomial) -> Polynomial:
        """
        Находит НОД двух многочленов на основе алгоритма Евклида.
        """
        # Копии многочленов
        pln1 = Polynomial({deg: coeff for deg, coeff in poly1.terms.items()})
        pln2 = Polynomial({deg: coeff for deg, coeff in poly2.terms.items()})

        while pln2.terms:  # Пока второй многочлен не пуст
            remainder = PolynomialOperations.MOD_PP_P(pln1, pln2)
            pln1, pln2 = pln2, remainder

        # Если старший коэффициент отрицательный, домножаем на -1
        lead_coeff = PolynomialOperations.LED_P_Q(pln1)
        if int(lead_coeff.numerator) < 0:
            pln1 = PolynomialOperations.MUL_PQ_P(pln1, Rational(Integer("-1"), Natural("1")))

        return pln1

    @staticmethod
    def DER_P_P(poly: Polynomial) -> Polynomial:
        """
        Находит производную многочлена.
        """
        derivative = Polynomial()

        for degree, coefficient in poly.terms.items():
            if degree > 0:  # Производная от константы равна 0
                new_coeff = RationalOperations.MUL_QQ_Q(coefficient, Rational(Integer(str(degree)), Natural("1")))
                derivative.terms[degree - 1] = new_coeff

        return derivative

    @staticmethod
    def NMR_P_P(poly: Polynomial) -> Polynomial:
        """
        Преобразует многочлен, устраняя кратные корни.
        """
        # Находим производную многочлена
        derivative = PolynomialOperations.DER_P_P(poly)

        # НОД между многочленом и его производной
        gcd_poly = PolynomialOperations.GCF_PP_P(poly, derivative)

        # Делим исходный многочлен на НОД
        result_poly = PolynomialOperations.DIV_PP_P(poly, gcd_poly)

        return result_poly

