/* nbody_kern.cl */
/* http://www.browndeertechnology.com/docs/BDT_OpenCL_Tutorial_NBody-rev3.html*/

__kernel void nbody(
   float dt1, float eps,
   __global float4* pos_old,
   __global float4* pos_new,
   __global float4* vel,
   __global float4* eterm,  /* rep, attract, lower, upper */
   __local float4* pblock
   __local float4* eblock
) {
   const float4 dt = (float4)(dt1,dt1,dt1,0.0f);

   int gti = get_global_id(0);
   int ti = get_local_id(0);
   
   int n = get_global_size(0);
   int nt = get_local_size(0);
   int nb = n/nt;

   int eoff = gti*n

   float4 p = pos_old[gti];
   float4 v = vel[gti];

   float4 a = (float4)(0.0f,0.0f,0.0f,0.0f);

    for(int jb=0; jb < nb; jb++) { /* Foreach block ... */

       pblock[ti] = pos_old[jb*nt+ti]; /* Cache x, y, z, w other particle */
       eblock[ti] = eterm[eoff+jb*nt+ti]; /* Cache energy terms */
       
       barrier(CLK_LOCAL_MEM_FENCE); /* Wait for others in the work-group */

       for(int j=0; j<nt; j++) { /* For CACHED particle positions ... */
          float4 p2 = pblock[j]; /* Read a cached particle position */
          float4 e = eblock[j]
          float4 d = p2 - p;
          
          /*  if bound do dist restraint, otherwise repulsive */
          /*  array of connectivity strength */
          /*  array of restraint bounds */
          /*  waterbath the velocities */
          
          float r2 = 
          
          
          float invr = rsqrt(d.x*d.x + d.y*d.y + d.z*d.z + eps);
          float f = p2.w*invr*invr*invr;
          
          a += f*d; /* Accumulate acceleration */
       }

       barrier(CLK_LOCAL_MEM_FENCE); /* Wait for others in work-group */
    }


    p += dt*v + 0.5f*dt*dt*a;
    v += dt*a;

    pos_new[gti] = p;
    vel[gti] = v;

 }

