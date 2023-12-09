'''
Any formatting details that could be used when making plots
'''


#########
# Colorbar number formatting functions
#########

def reg1_fmt(x, pos):
    return f'{x:.1f}'

def reg3_fmt(x, pos):
    return f'{x:.3f}'

def reg5_fmt(x, pos):
    return f'{x:.5f}'

def scinot2_fmt(x, pos):
    '''
    https://stackoverflow.com/questions/25983218/scientific-notation-colorbar
    '''
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def scinot3_fmt(x, pos):
    '''
    https://stackoverflow.com/questions/25983218/scientific-notation-colorbar
    '''
    a, b = '{:.3e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)