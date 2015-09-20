def kramer_solve_SLE_2(a11, a12, a21, a22, b1, b2):
    """Solves System of Linear Equations rank==2 using Kramer's rule"""
    _delta = a11 * a22 - a12 * a21
    _delta1 = b1 * a22 - b2 * a12
    _delta2 = b2 * a11 - b1 * a21
    
    x = _delta1 / _delta
    y = _delta2 / _delta
    return (x, y)

