'''
This module is a collection of functions for working with vectors.
'''
from math import sqrt, acos, fabs, modf, pi
from decimal import Decimal, getcontext

# Decimal precision
getcontext().prec = 30


class Vector(object):
    '''
    This class is a collection of functions for working with vectors.
    '''
    # Error messages
    CANNOT_NORMALIZE_ZERO_VECTOR_MSG = 'Cannot normalize a zero vector.'
    NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG = 'No unique orthogonal component: \
                                        At least one of the vectors is a zero vector.'
    NO_UNIQUE_PARALLEL_COMPONENT_MSG = 'No unique parallel component: \
                                        At least one of the vectors is a zero vector.'
    ONLY_DEFINED_IN_TWO_THREE_DIMS_MSG = 'Cross product is only defined for R^2 and R^3.'

    def __init__(self, coords):
        '''
        Init new Vector with an Array of coordinates
        '''
        self.count = 0
        try:
            if not coords:
                raise ValueError
            self.coords = tuple([Decimal(x) for x in coords])
            self.dimension = len(coords)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')

    def __iter__(self):
        self.current = 0
        return self

    def __next__(self):
        if self.current >= len(self.coords):
            raise StopIteration
        else:
            current_value = self.coords[self.current]
            self.current += 1
            return current_value

    def __len__(self):
        return len(self.coords)

    def __getitem__(self, i):
        return self.coords[i]

    def __str__(self):
        return 'Vector: {}'.format(self.coords)

    def __eq__(self, v):
        return self.coords == v.coords

    def add(self, v):
        '''
        Add two vectors.
        
        :param Vector v: second vector
        :return: result of addition
        :rtype: Vector
        '''
        new_coords = [x + y for x, y in zip(self.coords, v.coords)]
        return Vector(new_coords)

    def substract(self, v):
        '''
        Substract two vectors.
        
        :param Vector v: second vector
        :return: results of substraction
        :rtype: Vector
        '''
        new_coords = [x - y for x, y in zip(self.coords, v.coords)]
        return Vector(new_coords)

    def scalar_multiply(self, c):
        '''
        Multiply by a scalar.
        
        :param float c: scalar
        :return: scalar product
        :rtype: Vector
        '''
        new_coords = [x * Decimal(c) for x in self.coords]
        return Vector(new_coords)

    def magnitude(self):
        '''
        Find vector's magnitude.

        :return: magnitude of the vector
        :rtype: Decimal
        '''
        return Decimal(sqrt(sum([x**2 for x in self.coords])))

    def normalize(self):
        '''
        Normalize vector (aka find unit vector).
        
        :return: unit vector
        :rtype: Vector
        :raises Exception: if the vectors is a zero vector
        '''
        try:
            return self.scalar_multiply((Decimal('1.0') / self.magnitude()))
        except ZeroDivisionError:
            raise Exception(self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG)

    def dot_multiply(self, v):
        '''
        Dot multiply two vectors.
        
        :param Vector v: second vector
        :return: dot product
        :rtype: Decimal
        '''
        return Decimal(sum([x * y for x, y in zip(self.coords, v.coords)]))

    def angle(self, v, unit='rad'):
        '''
        Find angle between two vectors.

        :param Vector v: second vector
        :param str unit: 'rad' (radians) or 'deg' (degrees)
        :param float tolerance: precision tolerance
        :return: angle in radians or degrees
        :rtype: Decimal
        :raises Exception: if one of the vectors is a zero vector
        '''
        try:
            self_norm = self.normalize()
            v_norm = v.normalize()
            cosang = min(1, max(self_norm.dot_multiply(v_norm), -1))
            angle = acos(cosang)
            return {
                'rad': angle,
                'deg': angle * 180 / pi
            }.get(unit, angle)

        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception('Cannot compute an angle with a zero vector')
            else:
                raise e

    def is_parallel(self, v):
        '''
        Check if two vectors are parallel.
        
        :param Vector v: second vector
        :return: whether the vectors are parallel
        :rtype: Boolean
        '''
        return (self.is_zero()
                or v.is_zero()
                or self.angle(v) == 0
                or self.angle(v) == pi)

    def is_orthogonal(self, v, tolerance=1e-10):
        '''
        Check if two vectors are orthogonal
        
        :param Vector v: second vector
        :param float tolerance: precision tolerance
        :return: whether the vectors are orthogonal
        :rtype: Boolean
        '''
        return abs(self.dot_multiply(v)) < tolerance

    def is_zero(self, tolerance=1e-10):
        '''
        Check if zero vector.
        
        :param float tolerance: precision tolerance
        :return: whether the vector is a zero vector
        :rtype: Boolean
        '''
        return self.magnitude() < tolerance

    def project(self, basis):
        '''
        Find parallel component (aka, project onto basis vector).
        
        :param Vector basis: basis vector
        :param float tolerance: precision tolerance
        :return: parallel component
        :rtype: Vector
        :raises Exception: if one of the vectors is a zero vector
        '''
        try:
            u = basis.normalize()
            weight = self.dot_multiply(u)
            return u.scalar_multiply(weight)

        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MSG)
            else:
                raise e

    def orthogonal_component(self, basis):
        '''
        Find orthogonal component.
        
        :param Vector basis: basis vector
        :param float tolerance: precision tolerance
        :return: orthogonal component
        :rtype: Vector
        :raises Exception: if one of the vectors is a zero vector
        '''
        try:
            projection = self.project(basis)
            return self.substract(projection)

        except Exception as e:
            if str(e) == self.NO_UNIQUE_PARALLEL_COMPONENT_MSG:
                raise Exception(self.NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG)
            else:
                raise e

    def cross_multiply(self, v):
        '''
        Find cross product (vector orthogonal to both vectors)

        :param Vector v: second vector
        :return: cross product
        :rtype: Vector
        :raises Exception: if vectors not in R^2 or R^3
        '''
        if self.dimension > 3 or self.dimension < 2 or v.dimension > 3 or v.dimension < 2:
            raise Exception(self.ONLY_DEFINED_IN_TWO_THREE_DIMS_MSG)
        elif self.dimension == 2 or v.dimension == 2:
            self_embedded_in_R3 = Vector(self.coords + ('0',))
            v_embedded_in_R3 = Vector(v.coords + ('0',))
            return self_embedded_in_R3.cross_multiply(v_embedded_in_R3)

        # a_2 * b_3 - b_2 * a_3
        cross_x = (self.coords[1] * v.coords[2]) - \
            (v.coords[1] * self.coords[2])
        # a_3 * b_1 - b_3 * a_1
        cross_y = (self.coords[2] * v.coords[0]) - \
            (v.coords[2] * self.coords[0])
        # a_1 * b_2 - b_1 * a_2
        cross_z = (self.coords[0] * v.coords[1]) - \
            (v.coords[0] * self.coords[1])

        return Vector([cross_x, cross_y, cross_z])

    def parallelogram_area(self, v):
        '''
        Calculates area of parallelogram spanned by 3D vectors.

        :param Vector v: second vector
        :return: parallelogram area
        :rtype: Decimal
        :raises Exception: if vectors not in R^2 or R^3
        '''
        if self.dimension > 3 or self.dimension < 2 or v.dimension > 3 or v.dimension < 2:
            raise Exception(self.ONLY_DEFINED_IN_TWO_THREE_DIMS_MSG)

        return self.cross_multiply(v).magnitude()

    def triangle_area(self, v):
        '''
        Calculates area of triangle half of a parallelogram spanned by 3D vectors.

        :param Vector v: second vector
        :return: triangle area
        :rtype: Decimal
        :raises Exception: if vectors not in R^2 or R^3
        '''
        if self.dimension > 3 or self.dimension < 2 or v.dimension > 3 or v.dimension < 2:
            raise Exception(self.ONLY_DEFINED_IN_TWO_THREE_DIMS_MSG)

        return self.parallelogram_area(v) / 2
