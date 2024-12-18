shanks(n, r, s) = {
  smodn = Mod(s, n);
  rmodn = Mod(r, n);
  my(
    wk = Mod(0, n),
    wkm1 = Mod(3, n),
    wkp1 = Mod(2*r, n),
    wnegk = -rmodn / smodn,
    wnegkm1 = (rmodn/smodn)^2,
    wnegkp1= Mod(3, n),
    spowk = smodn,
    twosmodn=2*smodn,
    invsmodn=Mod(1,n)/smodn,
    twoinvsmodn=2*invsmodn,
    spownegk=Mod(1,n)/spowk;
  ); 

  my(bits = binary(n), k = 1); 

  for (j = 2, #bits,
    
    my(w2k = wk^2 - 2 * spowk * wnegk);    
                
    my(w2kplus2 = wkp1^2 - spowk*twosmodn * wnegkm1);  
       
    my(w2kminus2 = wkm1^2 - spowk * twoinvsmodn * wnegkp1);  
        
    my(wneg2k = wnegk^2 - 2 * spownegk * wk);    
          
    my(wneg2kplus2 = wnegkp1^2 - spownegk*twosmodn * wkm1);    
    
    my(wneg2kminus2 = wnegkm1^2 - (spownegk*twoinvsmodn) * wkp1);
    


        my(w2kminus1 = (w2kplus2 - rmodn * w2k) * invsmodn);    
           
        my(w2kplus1 = rmodn * w2kminus1 + smodn * w2kminus2);   

        my(wneg2kminus1 = (wneg2kplus2 - rmodn * wneg2k) * invsmodn);  
        my(wneg2kplus1 = rmodn * wneg2kminus1 + smodn * wneg2kminus2);





    if (bits[j], 
      
       
  
        
        wkm1 = w2k;
        wk = w2kplus1;
        wkp1 = w2kplus2;
        spowk=smodn*spowk^2;
        spownegk=(spownegk^2)/smodn;

        wnegkm1 = wneg2kminus2;
        wnegk = wneg2kminus1;
        wnegkp1 = wneg2k;
        k=2*k+1;
      );


    if (!bits[j],
      
       
        wkm1 = w2kminus1;
        wk = w2k;
        wkp1 = w2kplus1;
        spowk=spowk^2;
        spownegk=spownegk^2;

        wnegkm1 = wneg2kminus1;
        wnegk = wneg2k;
        wnegkp1 = wneg2kplus1;
        k=2*k;
      );
    
      );

    answer=matrix(7,1);
    answer[1,1]=wkp1;
    answer[2,1]=wk;
    answer[3,1]=wkm1;
    answer[4,1]=spowk;
    answer[5,1]=wnegk;
    answer[6,1]=wnegkm1;
    answer[7,1]=wnegkp1;
    return(answer);

  
};










isfrob3select(n) =
{
    if((n==1), return(0));
    if(((n%2)==0)&&(n>2), return(0));
    if(ispower(n), return(0));

    my(conductors = [7, 13, 19, 37, 61, 67, 73, 79, 97, 103, 139, 151, 163, 181, 193, 199, 211, 241, 271, 313, 331, 337, 349, 367, 373, 379, 409, 421, 463, 487, 523, 541, 547, 571, 577, 607, 613, 619, 631, 661, 673, 709, 751, 757, 769, 787, 823, 829, 853, 859, 877, 883, 907, 937, 967, 991, 1009, 1033, 1039, 1063, 1087, 1117, 1123, 1129, 1153, 1171, 1201, 1213, 1231, 1237, 1249, 1279, 1291, 1297, 1303, 1321, 1381, 1429, 1447, 1453, 1483, 1489, 1531, 1543, 1549, 1567, 1609, 1621, 1663, 1669, 1693, 1741, 1747, 1759, 1783, 1861, 1867, 1873, 1879, 1951],

       s_values = [7, 13, 19, 37, 183, 201, 219, 79, 97, 309, 139, 453, 163, 905, 579, 995, 1055, 1205, 813, 313, 2317, 2359, 349, 1101, 2611, 1895, 2045, 2947, 3241, 3409, 1569, 3787, 4923, 3997, 5193, 607, 1839, 5571, 3155, 1983, 4711, 709, 5257, 6813, 3845, 7083, 9053, 9119, 7677, 9449, 877, 6181, 9977, 937, 8703, 2973, 9081, 7231, 5195, 1063, 7609, 3351, 12353, 1129, 14989, 15223, 8407, 15769, 16003, 13607, 11241, 14069, 6455, 16861, 11727, 3963, 17953, 7145, 18811, 10171, 19279, 1489, 22965, 4629, 23235, 1567, 24135, 4863, 11641, 15021, 22009, 22633, 19217, 26385, 5349, 27915, 5601, 20603, 16911, 33167]);

    for (i = 1, #conductors, 
        my(N = conductors[i], r = N, s = s_values[i], rootD=sqrtint(4*r^3-27*s^2));  
        if((gcd(n,N)>1)&& (n>r), return(0)); 
        if((gcd(n,s)>1)&& (n>gcd(n,s)), return(0));   
        my(witness = n % N);
        while (!isprime(witness), witness += N);
        
        if (polisirreducible(Mod(1, witness)*X^3 - Mod(r, witness)*X - Mod(s, witness)),           
            my(M = matrix(3, 3), v3 = matrix(7, 1), v=matrix(3,1));
            M[1, 1] = Mod(0, n); M[1, 2] = Mod(r, n); M[1, 3] = Mod(s, n);
            M[2, 1] = Mod(1, n); M[2, 2] = M[1, 1]; M[2, 3] = M[1, 1];
            M[3, 1] = M[1, 1]; M[3, 2] = Mod(1, n); M[3, 3] = M[1, 1];
            v3= shanks(n, r, s); 
            v[1, 1] = v3[1, 1];
            v[2, 1] = v3[2, 1];
            v[3, 1] = v3[3, 1];        
            my(output = M * v);           
            my(test1 = (output[2, 1] == Mod(-r, n)),
               test2 = (output[3, 1] == Mod(0, n)),
               test4 = (v3[4,1] == Mod(s, n)),
               test5 = (v3[5,1] == Mod(-r,n)/Mod(s,n)),
               test6 = (v3[6,1] == Mod(0,n)),
               


               a=((-3*s + rootD)/2)%n, b=((-3*s - rootD)/2)%n, 
               test7=((v3[7,1]==Mod(a,n)/Mod(s,n))||(v3[7,1]==Mod(b,n)/Mod(s,n))),
               test3=((output[1,1]==Mod(a,n))||(output[1,1]==Mod(b,n))));
            if (test1 && test2 && test3 && test4&&test5&&test6&&test7, 

            
            
            
            
              
               
               return(1); 
               );
               if(!(test1 && test2 && test3 && test4&&test5&&test6&&test7),
                   return(0);
                 )
             
        )
    );

    return(-1); 
}
