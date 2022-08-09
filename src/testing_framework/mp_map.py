from multiprocessing import Pool
import time

class a:
    def __init__(self, in1, in2, in3): #This is like the named tuple
        self.num1 = in1
        self.num2 = in2
        self.num3 = in3
    
    def __call__(self, o): #This is like id_spectrum
        time.sleep(1)
        return o + self.num1

if __name__ == '__main__':
    x_tot = [a(4,5,6), a(12,13,14), a(21,22,23)]
    y_tot = []
    
    p = Pool(1)
    now = time.time()
    z = [i for i in range(10)]
    x = a(4,5,6)
    print(x(5))
    y = p.map(x_tot, z) #function can only take 1 input so make object
    y_tot.append(y)
    now = time.time() - now
    print(now)
    print(y_tot)