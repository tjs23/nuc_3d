/* nbody_kern.cl */
/* http://www.browndeertechnology.com/docs/BDT_OpenCL_Tutorial_NBody-rev3.html*/


/* Restrains stored as flat array of partners, with corresponding index, count array */

__kernel void rpd(float dt1, float rest_const, float rep_const, float rep_dist_2, float eps,
                  __global float4* pos_old, __global float4* pos_new,
                  __global float4* vel,     __global int* rest_pairs,
                  __global int2* rest_ic,   __global float2* rest_dists,
                  __local float4* pblock)
{
                  
   const float4 dt = (float4)(dt1,dt1,dt1,0.0f);

   int i_g = get_global_id(0);
   int i_l = get_local_id(0);

   int n = get_global_size(0);
   int n_l = get_local_size(0);
   int n_b = n/n_l;  
   
   float4 p = pos_old[i_g];
   float4 v = vel[i_g];
   float4 a = (float4)(0.0f,0.0f,0.0f,0.0f);
   int2 index_count = rest_ic[i_g]
   
   /* Restraints for this particle */
   
   for (int j=0; j<index_count.y; j++) {
     int k = rest_pairs[index_count.x+j]   /* Index of other particle */
     
     float2 rd = rest_dists[k];  /* Min and max restraint dists */
     float4 p2 = pos_old[k];     /* Other coords */
     float4 dp = p2 - p;         /* Delta coords */
     
     float r = sqrt(dp.x*dp.x + dp.y*dp.y + dp.z*dp.z);
     float f = 2 * rest_const * (max(rd.x - r, 0.0f) - max(r - rd.y, 0.0f)) /* Flat bottomed harmonic */
     
     a += f*dp; /* Accumulate acceleration */
   
   }
   
   /* All vs all repulsions */
   
   for (int block=0; block < n_b; block++) { /* For each block ... */

      pblock[i_l] = pos_old[block * n_l + i_l]; /* Cache ONE particle position */
      barrier(CLK_LOCAL_MEM_FENCE);             /* Wait for others in the work-group */

      for(int j=0; j<n_l; j++) { /* For ALL cached particle positions ... */
         float4 p2 = pblock[j];  /* Read a cached particle position,  p2.w is mass */
         float4 dp = p2 - p;
         
         float r2 = dp.x*dp.x + dp.y*dp.y + dp.z*dp.z + eps;
         
         float f = 4 * rep_const * max(rep_dist_2 - r2, 0.0f);
         
         a += f*dp; /* Accumulate acceleration */
      }

      barrier(CLK_LOCAL_MEM_FENCE); /* Wait for others in work-group */
   }
  
   /* Add water bath to accel before updating motion */
   
   p += dt*v + 0.5f*dt*dt*a;
   v += dt*a;

   pos_new[i_g] = p;
   vel[i_g] = v;

   /* Update velocity for water bath */
 
   decay = exp(-Langevin_gamma * 0.5f*dt)
   R = Langevin_gamma * BOLTZMANN_K * tRef / 0.5f*dt
 
   a = f/mass + sqrt(R/mass) * normSampl
 
 
   v = (decay * v) + (0.5f * dt * a)
   p += dt*v
  
}

