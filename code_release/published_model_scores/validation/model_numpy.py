"""CRISPR-Net scoring forward pass; NumPy arithmetic, h5py weight reader.

Architecture: upstream commit 0843c485b1d476f899bedb93a7dd2df42547c1e1.
All arithmetic is float64; stored float32 parameters are promoted without fitting.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '4'
import numpy as np
import h5py


def encode(pairs):
    """Published seven-channel OR/direction encoder, with separate left padding."""
    code = {'A': (1,0,0,0,0), 'T': (0,1,0,0,0), 'G': (0,0,1,0,0),
            'C': (0,0,0,1,0), '_': (0,0,0,0,1), '-': (0,0,0,0,0)}
    direction = {'A':5, 'G':4, 'C':3, 'T':2, '_':1}
    result = np.zeros((len(pairs),24,7), dtype=np.float64)
    for n,(on,off) in enumerate(pairs):
        if not (0 < len(on) <= 24 and 0 < len(off) <= 24):
            raise ValueError('Sequences must contain 1..24 characters')
        on,off = on.rjust(24,'-'),off.rjust(24,'-')
        for t,(a,b) in enumerate(zip(on,off)):
            if a == 'N': a = b
            result[n,t,:5] = np.bitwise_or(code[a],code[b])
            if a != '-' and b != '-' and direction[a] != direction[b]:
                result[n,t,5 if direction[a] > direction[b] else 6] = 1
    return result


class NumpyModel:
    def __init__(self, path, dtype=np.float64):
        self.w = {}
        with h5py.File(path,'r') as f:
            f.visititems(lambda k,v: self.w.__setitem__(k, v[()].astype(dtype))
                         if isinstance(v,h5py.Dataset) else None)
        shapes = {f'conv2d_{i}/conv2d_{i}/kernel:0': (1,k,7,10)
                  for i,k in enumerate((1,2,3,5),1)}
        shapes.update({'dense_1/dense_1/kernel:0':(720,80),
                       'dense_2/dense_2/kernel:0':(80,20),
                       'main_output/main_output/kernel:0':(20,1)})
        for key,shape in shapes.items():
            if self.w[key].shape != shape: raise ValueError((key,self.w[key].shape))

    def lstm(self,x,direction):
        prefix = f'bidirectional_1/bidirectional_1/{direction}_LSTM_out/'
        w,r,b = (self.w[prefix+k+':0'] for k in ('kernel','recurrent_kernel','bias'))
        n,t,_ = x.shape
        h,c = np.zeros((n,15),dtype=x.dtype),np.zeros((n,15),dtype=x.dtype)
        z = x @ w + b
        out = np.empty((n,t,15),dtype=x.dtype)
        order = range(t) if direction == 'forward' else range(t-1,-1,-1)
        for j in order:
            gates = z[:,j] + h @ r
            i,f,g,o = np.split(gates,4,axis=1)
            i = np.clip(.2*i+.5,0,1)
            f = np.clip(.2*f+.5,0,1)
            o = np.clip(.2*o+.5,0,1)
            c = f*c + i*np.tanh(g)
            h = o*np.tanh(c)
            out[:,j] = h
        return out

    def predict(self,x):
        if x.shape[1:] != (24,7): raise ValueError(x.shape)
        branches = [x]
        for i in range(1,5):
            prefix = f'conv2d_{i}/conv2d_{i}/'
            w,b = self.w[prefix+'kernel:0'][0],self.w[prefix+'bias:0']
            k = len(w)
            left = (k-1)//2
            padded = np.pad(x,((0,0),(left,k-1-left),(0,0)))
            y = np.broadcast_to(b,(len(x),24,10)).copy()
            for j in range(k): y += padded[:,j:j+24] @ w[j]
            branches.append(np.maximum(y,0))
        mixed = np.concatenate(branches,axis=-1)
        y = np.concatenate([self.lstm(mixed,'forward'),
                            self.lstm(mixed,'backward')],axis=-1).reshape(len(x),720)
        for name in ('dense_1','dense_2'):
            p = f'{name}/{name}/'
            y = np.maximum(y @ self.w[p+'kernel:0'] + self.w[p+'bias:0'],0)
        logits = (y @ self.w['main_output/main_output/kernel:0'] +
                  self.w['main_output/main_output/bias:0']).ravel()
        # Stable sigmoid, retaining small values lost by the historical float32 backend.
        ans = np.empty_like(logits)
        pos = logits >= 0
        ans[pos] = 1/(1+np.exp(-logits[pos]))
        e = np.exp(logits[~pos]); ans[~pos] = e/(1+e)
        return ans


def batches(model,pairs,size=512):
    return np.concatenate([model.predict(encode(pairs[j:j+size]))
                           for j in range(0,len(pairs),size)])
