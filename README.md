# utl-mle-symbolic-solution-for-mu-and-sigma-of-normal-pdf-using-sympy
Symbolic solution for population  mu and sigma  of normal pdf using sympy 
    %let pgm=utl-mle-symbolic-solution-for-mu-and-sigma-of-normal-pdf-using-sympy;

    Symbolic solution for population  mu and sigma  of normal pdf using sympy

    gihub
    https://tinyurl.com/kxhnk4sv
    https://github.com/rogerjdeangelis/utl-mle-symbolic-solution-for-mu-and-sigma-of-normal-pdf-using-sympy

    Two Solutions

          1 python sympy (proof by induction?)
          2 theory

     We will show that n independent normal random variables from a population of size n that

       mu    = sum(x1,x2,..xn)/n

                     ________________
                    /     /1       2
                   / sum | (xi-mu)
       sigma =    /       \n           (sometimes called the unadjusted standard deviation)
                \/  -----------------
                          n
    /*                   _
    (_)_ __  _ __  _   _| |_
    | | `_ \| `_ \| | | | __|
    | | | | | |_) | |_| | |_
    |_|_| |_| .__/ \__,_|\__|
            |_|
    */

    /*****************************************************************************************************************************************/
    /*                                    |                                                               |                                  */
    /*              INPUT                 |                 PROCESS                                       |       OUTPUT                     */
    /*                                    |                                                               | Proof by Induction using sympy   */
    /* X1, X2, ..., Xn be a random sample | %utl_pybegin;                                                 |                                  */
    /* from a normal distribution with    | parmcards4;                                                   | Sample size 3 only to show       */
    /* unknown mean mu and unknown        | from sympy import symbols, exp, pi, Integral \                | process & prove symbolic result  */
    /* variance sigma                     |  , solve, diff, log, sqrt, simplify, factor, pprint           | symbolic result                  */
    /*                                    | mu, sigma, x = symbols('mu sigma x', positive=True)           |                                  */
    /*                                    | n = symbols('n', positive=True, integer=True)                 | From Sympy estimates of mu       */
    /*                               2    | x1,x2,x3 = symbols('x1 x2 x3',positive=True)                  |                                  */
    /*                     -(-mu + x)     | norm_pdf=1/(sqrt(2*pi*sigma**2))*exp(-(x-mu)**2/(2*sigma**2)) | mu=(x1+x2+x3)/3                  */
    /*                     ------------   | # pprint(norm_pdf)                                            |                                  */
    /*                              2     | log_likelihood=sum(log(norm_pdf.subs(x,x_i))                  | sigma=sqrt(3)*                   */
    /*                ___    2*sigma      |    for x_i in(x1,x2,x3))                                      | sqrt(3*mu**2                     */
    /*              \/ 2 *e               | #pprint(log_likelihood);                                      |   -2*mu*x1                       */
    /* normalpdf =  -------------------   | derivative = diff(log_likelihood, mu)                         |   -2*mu*x2                       */
    /*                     ____           | mu_mle_equation = derivative.simplify()                       |   -2*mu*x3                       */
    /*                 2*\/ pi *sigma     | #pprint(derivative)                                           |   +x1**2                         */
    /*                                    | mu_mle = solve(mu_mle_equation, mu)[0]                        |   +x2**2                         */
    /*  Verify with                       | print(f"The MLE of mu is: {mu_mle}")                          |   +x3**2)/3;                     */
    /*                                    | derivative = diff(log_likelihood, sigma)                      |                                  */
    /*    x1=10                           | sig_mle_equation = derivative.simplify()                      |  Which simplifies to             */
    /*    x2=20                           | #pprint(sig_mle_equation)                                     |  see below                       */
    /*    x3=10                           | sig_mle = solve(sig_mle_equation, sigma)[0]                   |                                  */
    /*                                    | pprint(simplify(sig_mle))                                     |              _____________       */
    /*  data nums;                        | print(simplify(sig_mle))                                      |            /    /1      2        */
    /*  input x @@;                       | #pprint(f"The MLE of sigma is: {sig_mle}")                    |           /sum | (xi-mu)         */
    /*  cards4;                           | ;;;;                                                          |  sigma=  /       \3              */
    /*  10 20 10                          | %utl_pyend;                                                   |        \/  ---------------       */
    /*  ;;;;                              |                                                               |                   3              */
    /*  run;quit;                         | Lets test using known formulas                                |                                  */
    /*                                    |    x1=10                                                      |  Check                           */
    /*                                    |    x2=20                                                      |  x1=10                           */
    /*                                    |    x3=10                                                      |  x2=20                           */
    /*                                    | data nums;                                                    |  x3=10                           */
    /*                                    | input x @@;                                                   |                                  */
    /*                                    | cards4;                                                       |   Proc means                     */
    /*                                    | 10 20 10                                                      |                                  */
    /*                                    | ;;;;                                                          |          Mean    Std Dev         */
    /*                                    | run;quit;                                                     |    ---------------------         */
    /*                                    | * vardef=n gives us population stats;                         |    13.3333333  4.7140452         */
    /*                                    | proc means data=nums vardef=N mean stddev;                    |    ---------------------         */
    /*                                    | run;quit;                                                     |                                  */
    /*                                    |                                                               |  Sympy                           */
    /*                                    |       Mean         Std Dev                                    |                                  */
    /*                                    | --------------------------                                    |  AGREES WITH PROC MEANS          */
    /*                                    | 13.3333333       4.7140452                                    |  ----------------------          */
    /*                                    | --------------------------                                    |  pop sigma: 4.7140452079         */
    /*                                    |                                                               |  mean: 13.333333333              */
    /*                                    | Now lets use MLE results                                      |                                  */
    /*                                    |                                                               |                                  */
    /*                                    | data chk;                                                     |                                  */
    /*                                    | array xs[3] x1-x3                                             |                                  */
    /*                                    |  (10 20 10);                                                  |                                  */
    /*                                    | mu=mean(of xs[*]);                                            |                                  */
    /*                                    | put "mean: " mu;                                              |                                  */
    /*                                    | sigma=sqrt(3)*                                                |                                  */
    /*                                    |   sqrt(3*mu**2                                                |                                  */
    /*                                    |     -2*mu*x1                                                  |                                  */
    /*                                    |     -2*mu*x2                                                  |                                  */
    /*                                    |     -2*mu*x3                                                  |                                  */
    /*                                    |     +x1**2                                                    |                                  */
    /*                                    |     +x2**2                                                    |                                  */
    /*                                    |     +x3**2)/3;                                                |                                  */
    /*                                    | put "pop sigma: " sigma;                                      |                                  */
    /*                                    | run;quit;                                                     |                                  */
    /*                                    |                                                               |                                  */
    /*                                    | AGREES WITH PROC MEANS                                        |                                  */
    /*                                    | ----------------------                                        |                                  */
    /*                                    | mean: 13.333333333                                            |                                  */
    /*                                    | pop sigma: 4.7140452079                                       |                                  */
    /*                                    |                                                               |                                  */
    /*                                    |                                                               |                                  */
    /*****************************************************************************************************************************************/

    /*                   _
    (_)_ __  _ __  _   _| |_
    | | `_ \| `_ \| | | | __|
    | | | | | |_) | |_| | |_
    |_|_| |_| .__/ \__,_|\__|
            |_|
    */

     X1, X2, ..., Xn be a random sample
     from a normal distribution with
     unknown mean mu and unknown
     variance sigma


                                   2
                         -(-mu + x)
                         ------------
                                  2
                    ___    2*sigma
                  \/ 2 *e
     normalpdf =  -------------------
                         ____
                     2*\/ pi *sigma

     Verify with

        x1=10
        x2=20
        x3=10

      data nums;
      input x @@;
      cards4;
      10 20 10
      ;;;;
      run;quit;

    /*
    / |  ___ _   _ _ __ ___  _ __  _   _   _ __  _ __ ___   ___ ___  ___ ___
    | | / __| | | | `_ ` _ \| `_ \| | | | | `_ \| `__/ _ \ / __/ _ \/ __/ __|
    | | \__ \ |_| | | | | | | |_) | |_| | | |_) | | | (_) | (_|  __/\__ \__ \
    |_| |___/\__, |_| |_| |_| .__/ \__, | | .__/|_|  \___/ \___\___||___/___/
             |___/          |_|    |___/  |_|
    */

    %utl_pybegin;
    parmcards4;
    from sympy import symbols, exp, pi, Integral \
     , solve, diff, log, sqrt, simplify, factor, pprint
    mu, sigma, x = symbols('mu sigma x', positive=True)
    n = symbols('n', positive=True, integer=True)
    x1,x2,x3 = symbols('x1 x2 x3',positive=True)
    norm_pdf=1/(sqrt(2*pi*sigma**2))*exp(-(x-mu)**2/(2*sigma**2))
    # pprint(norm_pdf)
    log_likelihood=sum(log(norm_pdf.subs(x,x_i))
       for x_i in(x1,x2,x3))
    #pprint(log_likelihood);
    derivative = diff(log_likelihood, mu)
    mu_mle_equation = derivative.simplify()
    #pprint(derivative)
    mu_mle = solve(mu_mle_equation, mu)[0]
    print(f"The MLE of mu is: {mu_mle}")
    derivative = diff(log_likelihood, sigma)
    sig_mle_equation = derivative.simplify()
    #pprint(sig_mle_equation)
    sig_mle = solve(sig_mle_equation, sigma)[0]
    #pprint(simplify(sig_mle))
    print(f"The MLE of sigma is: {sig_mle}")
    #pprint(f"The MLE of sigma is: {sig_mle}")
    ;;;;
    %utl_pyend;


    * CHECK SYMPY RESULTS;

    data nums;
    input x @@;
    cards4;
    10 20 10
    ;;;;
    run;quit;
    * vardef=n gives us population stats;
    proc means data=nums vardef=N mean stddev;
    run;quit;

    /*
          Mean         Std Dev
    --------------------------
    13.3333333       4.7140452
    --------------------------
    */

    Now lets use MLE results

    data chk;
    array xs[3] x1-x3
     (10 20 10);
    mu=mean(of xs[*]);
    put "mean: " mu;
    sigma=sqrt(3)*
      sqrt(3*mu**2
        -2*mu*x1
        -2*mu*x2
        -2*mu*x3
        +x1**2
        +x2**2
        +x3**2)/3;
    put "pop sigma: " sigma;
    run;quit;

    AGREES WITH PROC MEANS
    ----------------------
    mean: 13.333333333
    pop sigma: 4.7140452079

    /*           _               _
      ___  _   _| |_ _ __  _   _| |_
     / _ \| | | | __| `_ \| | | | __|
    | (_) | |_| | |_| |_) | |_| | |_
     \___/ \__,_|\__| .__/ \__,_|\__|
                    |_|
    */

    /**************************************************************************************************************************/
    /*                                                                                                                        */
    /*                                                                                                                        */
    /* The MLE of mu is: x1/3 + x2/3 + x3/3                                                                                   */
    /* The MLE of sigma is: sqrt(3)*sqrt(3*mu**2 - 2*mu*x1 - 2*mu*x2 - 2*mu*x3 + x1**2 + x2**2 + x3**2)/3                     */
    /*                                                                                                                        */
    /*  Simplifies to                                                                                                         */
    /*                                                                                                                        */
    /*      Note    2                2      2      2                                                                          */
    /*            x1  - 2*mu*x1  + mu   = (x1 - mu)                                                                           */
    /*              2                2      2      2                                                                          */
    /*            x2  - 2*mu*x2  + mu   = (x3 - mu)                                                                           */
    /*              2                2      2      2                                                                          */
    /*            x3  - 2*mu*x3  + mu   = (x1 - mu)                                                                           */
    /*                                                                                                                        */
    /*      Therfore                                                                                                          */
    /*                   _____________                                                                                        */
    /*                 /    /1      2                                                                                         */
    /*                /sum | (xi-mu)                                                                                          */
    /*       sigma=  /       \3                                                                                               */
    /*             \/  ---------------                                                                                        */
    /*                        3                                                                                               */
    /*                                                                                                                        */
    /**************************************************************************************************************************/

    /*___    _   _
    |___ \  | |_| |__   ___  ___  _ __ _   _
      __) | | __| `_ \ / _ \/ _ \| `__| | | |
     / __/  | |_| | | |  __/ (_) | |  | |_| |
    |_____|  \__|_| |_|\___|\___/|_|   \__, |
                                       |___/
    */

    Demonstrate that the MLE estimates for mu and sigma for the normal pdf
    for  random independent samples opoulation size n

      For our random sample of size 3

      mu = (x1 + x2 + x3 )/s

                  _____________
                /    /1      2
               /sum | (xi-mu)
      sigma=  /       \3
            \/  ---------------
                       3


      by induction

      mu = (x1 + x2 .. xn )/n

                  _____________
                /    /1      2
               /sum | (xi-mu)
      sigma=  /       \n
            \/  ---------------
                       n

    Let X1, X2, ..., Xn be a random sample from a normal distribution with
      unknown mean µ and unknown variance s^2.
    For documentation purposes I am only using n=3 or x1, x2 and x3

    The pdf of the normal distribution

                                   2
                         -(-mu + x)
                         ------------
                                  2
                    ___    2*sigma
                  \/ 2 *e
     normalpdf =  -------------------
                         ____
                     2*\/ pi *sigma

    The likelihood function is the joint pdf of the sample.

    Because of independence the joint probability
    is just the produce of univariate pdfs.



                                   2                        2                      2
                        -(-mu + x1)              -(-mu + x2)             -(-mu + x3)
                        ------------             ------------            ------------
                                 2                        2                       2
                   ___    2*sigma           ___    2*sigma          ___    2*sigma
                 \/ 2 *e                  \/ 2 *e                 \/ 2 *e
     joint pdf = -------------------  *   -------------------  *  -------------------
                        ____                     ____                    ____
                    2*\/ pi *sigma           2*\/ pi *sigma          2*\/ pi *sigma


     We need to find the values mu and sigma that mazimize the probability of the joint pdf.
     ie maximize the probability

    Easier to do in we take the log of the liklihood.


                            /                  2 \      /                  2 \      /                  2 \
                            |       -(-mu + x1)  |      |       -(-mu + x2)  |      |       -(-mu + x3)  |
                            |       -------------|      |       -------------|      |       -------------|
                            |                 2  |      |                 2  |      |                 2  |
                            |  ___     2*sigma   |      |  ___     2*sigma   |      |  ___     2*sigma   |
                            |\/ 2 *e             |      |\/ 2 *e             |      |\/ 2 *e             |
    ljpdf=log joint pdf =log|--------------------| + log|--------------------| + log|--------------------|
                            |       ____         |      |       ____         |      |       ____         |
                            \   2*\/ pi *sigma   /      \   2*\/ pi *sigma   /      \   2*\/ pi *sigma   /


    Lets solve for the vales of x1, x2 and x3 that maximixe the probablity
    by setting the  partial derivative with respect to mu equal to 0


                       dljpdf     2*mu - 2*x1   2*mu - 2*x2   2*mu - 2*x3
    partial derivative --     = - ----------- - ----------- - -----------  = 0
                       du                  2             2             2
                                    2*sigma       2*sigma       2*sigma

    Lets solve for mu

    Note we can ignore the denominators and the 2 factor

        -mu+x1 - mu+x2 - mu+ x3

       -3*mu + (x1 + x2 +x3) = 0

        mu =(x1+x2+x3)/3

    Lets solve for sigma by setting the partial derivative of the log joint pdf to 0

                                      2                                        2     2     2     2
                       dljpdf     3*mu  - 2*mu*x1 - 2*mu*x2 - 2*mu*x3 - 3*sigma  + x1  + x2  + x3
    partial derivative ----- =    ----------------------------------------------------------------  = 0
                       dsigma                                      3
                                                              sigma
             dljpdf
    solving  --    =  0
             du
                      _______________________________________________________
               ___   /     2                                   2     2     2
              / 3 *\/  3*mu  - 2*mu*x1 - 2*mu*x2 - 2*mu*x3 + x1  + x2  + x3
     sigma =  ---------------------------------------------------------------
                                            3


     Rearanging

         Note    2                2      2      2
               x1  - 2*mu*x1  + mu   = (x1 - mu)
                 2                2      2      2
               x2  - 2*mu*x2  + mu   = (x3 - mu)
                 2                2      2      2
               x3  - 2*mu*x3  + mu   = (x1 - mu)

         Therfore

                        _______________
                 ___   /    1        2
                / 3 *\/  sum  (xi-mu)
                            3
                -----------------------
                          3

          Multiply numerator and denominator sqrt(3)
          Population standard deviation

                 _______________
                /    1        2
               /  sum  (xi-mu)
              /      n
            \/  ---------------
                      3

    /*              _
      ___ _ __   __| |
     / _ \ `_ \ / _` |
    |  __/ | | | (_| |
     \___|_| |_|\__,_|

    */
