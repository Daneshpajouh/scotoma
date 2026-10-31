const assert = require('node:assert/strict');
const {inflateSync} = require('node:zlib');
const {createHash} = require('node:crypto');
const {Worker} = require('node:worker_threads');
const threads = Math.min(24, require('node:os').availableParallelism());

const encoded = `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`;

const raw = inflateSync(Buffer.from(encoded.replace(/\s/g,''), 'base64'));
assert.equal(createHash('sha256').update(raw).digest('hex'),
 '2ab76785c7602e5c9bddfbcb135e1b1f1d3326637de2bd1b7836c08a96ceff4e');

const N=19, M=8388607, size=3**N;
const shared=new SharedArrayBuffer(4*size);
const a=new Uint32Array(shared), seen=new Set();
let sum=0;

for(let off=0;off<raw.length;off+=5){
 const index=raw.readUInt32LE(off), b=raw[off+4];
 let p=index, stars=0, nonzero=0, walk=0, peak=0;
 for(let j=0;j<N;j++){
  const c=p%3; p=Math.floor(p/3);
  stars+=(c===2); nonzero+=(c!==0);
  walk+=(c===1?1:c===0?-1:0); peak=Math.max(peak,walk);
 }
 assert(p===0 && b>0 && stars<=8 &&
        nonzero-peak===8 && !seen.has(index));
 seen.add(index);
 a[index]=b*2**(8-stars);
 sum+=b;
}
assert.equal(seen.size,564);
assert.equal(sum,30348);
assert(sum*256<M);

function kernel(){
 const {parentPort,workerData:d}=require('node:worker_threads');
 const a=new Uint32Array(d.shared), M=8388607, size=3**19, split=3**9;

 function degrees(n){
  const v=new Uint8Array(3**n);
  for(let i=1;i<v.length;i++)
   v[i]=v[Math.floor(i/3)]+(i%3!==0);
  return v;
 }
 const low=degrees(9),high=degrees(10);
 function red(v){
  v=(v&M)+Math.floor(v/8388608);
  return v>=M?v-M:v;
 }

 parentPort.on('message',q=>{
  const total=(q.op==='F'||q.op==='I')?size/3:size;
  const start=Math.floor(d.id*total/d.threads);
  const end=Math.floor((d.id+1)*total/d.threads);

  if(q.op==='F'||q.op==='I'){
   const s=q.stride;
   let p=Math.floor(start/s)*3*s+start%s;
   let left=s-start%s,i=start;
   while(i<end){
    const len=Math.min(left,end-i),stop=p+len;
    for(;p<stop;p++){
     const x=a[p],y=a[p+s],z=a[p+2*s];
     if(q.op==='F'){
      a[p]=red(x+y+2*z);
      a[p+s]=red(x+M-y);
      a[p+2*s]=red(x+y+2*M-2*z);
     }else{
      a[p]=red(x+2*y+z);
      a[p+s]=red(x+2*M-2*y+z);
      a[p+2*s]=red(x+M-z);
     }
    }
    i+=len; p+=2*s; left=s;
   }
  }else if(q.op==='D'){
   for(let i=start;i<end;i++){
    const k=low[i%split]+high[Math.floor(i/split)];
    a[i]=red(a[i]*q.lambda[k]);
   }
  }else{
   let max=-1,arg=-1,min=M;
   for(let i=start;i<end;i++){
    if(a[i]>max){max=a[i];arg=i;}
    min=Math.min(min,a[i]);
   }
   parentPort.postMessage({max,arg,min});
   return;
  }
  parentPort.postMessage(null);
 });
}

(async()=>{
 const pool=Array.from({length:threads},(_,id)=>
  new Worker('('+kernel.toString()+')()',
   {eval:true,workerData:{shared,id,threads}}));

 async function phase(q){
  return Promise.all(pool.map(w=>new Promise((resolve,reject)=>{
   const done=x=>{w.off('error',fail);resolve(x);};
   const fail=e=>{w.off('message',done);reject(e);};
   w.once('message',done);
   w.once('error',fail);
   w.postMessage(q);
  })));
 }

 function C(n,k){
  if(k<0||k>n)return 0;
  let v=1;
  for(let j=1;j<=k;j++)v=v*(n-j+1)/j;
  return v;
 }

 const lambda=Array.from({length:20},(_,k)=>{
  let v=0;
  for(let j=0;j<=8;j++)
   for(let l=0;l<=j;l++)
    v+=(-1)**l*3**(j-l)*C(k,l)*C(19-k,j-l);
  return (((v%M+M)%M)*256)%M;
 });

 const start=Date.now();
 for(let j=0;j<19;j++)await phase({op:'F',stride:3**j});
 await phase({op:'D',lambda});
 for(let j=18;j>=0;j--)await phase({op:'I',stride:3**j});

 const rr=await phase({op:'max'});
 const best=rr.reduce((a,b)=>a.max>b.max?a:b);
 assert.equal(best.max,402688);
 assert.equal(Math.min(...rr.map(r=>r.min)),0);

 let index=best.arg,centre='';
 for(let j=0;j<19;j++){
  centre+=index%3;
  index=Math.floor(index/3);
 }
 assert.equal(centre,'1001010000001000000');
 assert(256*sum>19*best.max);

 console.log({
  orbits:size,
  max_integer_score:best.max,
  centre,
  mass:sum+'/1573',
  strict_margin_numerator:sum-19*1573,
  seconds:(Date.now()-start)/1000
 });
 await Promise.all(pool.map(w=>w.terminate()));
})().catch(e=>{console.error(e);process.exit(1)});
