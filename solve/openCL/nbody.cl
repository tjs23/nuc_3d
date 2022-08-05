/* nbody_kern.cl */
/* http://www.browndeertechnology.com/docs/BDT_OpenCL_Tutorial_NBody-rev3.html*/

__kernel void nbody(
   float dt1, float eps,
   __global float4* pos_old,
   __global float4* pos_new,
   __global float4* vel,
   __local float4* pblock
) {
   const float4 dt = (float4)(dt1,dt1,dt1,0.0f);

   int i_g = get_global_id(0);
   int i_l = get_local_id(0);

   int n = get_global_size(0);
   int n_l = get_local_size(0);
   int n_b = n/n_l;

   float4 p = pos_old[i_g];
   float4 v = vel[i_g];

   float4 a = (float4)(0.0f,0.0f,0.0f,0.0f);

    for(int block=0; block < n_b; block++) { /* For each block ... */

       pblock[i_l] = pos_old[block * n_l + i_l]; /* Cache ONE particle position */
       barrier(CLK_LOCAL_MEM_FENCE);             /* Wait for others in the work-group */

       for(int j=0; j<n_l; j++) { /* For ALL cached particle positions ... */
          float4 p2 = pblock[j];  /* Read a cached particle position */
          float4 d  = p2 - p;
          float invr = rsqrt(d.x*d.x + d.y*d.y + d.z*d.z + eps);
          float f    = p2.w*invr*invr*invr;
          a += f*d; /* Accumulate acceleration */
       }

       barrier(CLK_LOCAL_MEM_FENCE); /* Wait for others in work-group */
    }


    p += dt*v + 0.5f*dt*dt*a;
    v += dt*a;

    pos_new[i_g] = p;
    vel[i_g] = v;

 }

