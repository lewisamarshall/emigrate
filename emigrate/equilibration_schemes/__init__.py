from Equilibrator import Equilibrator
from Fixed import Fixed
from VariablepH import VariablepH

equilibrators = {None: Equilibrator,
                 'fixed': Fixed,
                 'pH': VariablepH,
                 }

if __name__ == '__main__':
    print equilibrators
