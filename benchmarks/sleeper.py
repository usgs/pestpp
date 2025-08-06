import numpy as np 
import time

def take_nap():
    minutes = np.random.uniform(0.1,5)
    start = time.time()
    print(f'sleeping for {minutes:0.2f} minutes')
    time.sleep(minutes*60)
    print(f"wake up! It's been {(time.time()-start)/60} minutes")

if __name__=="__main__":
    take_nap()
