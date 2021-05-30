/* =============================================================================  
** Color graphing system - mike3 - 20.11.10 04:07
** Hi.
** I thought I'd post the code I use to generate the color graphs from Pari/GP.
** Here it is.
**
** Note: the output is in .PPM format. 
** You'll need something else to convert that to .PNG. (I use GIMP.)
** 
** Also, I might warn you: it takes a LONG time to graph a complex function
** with a significantly complicated calculation procedure, as it must be
** evaluated at every pixel of the graph.
** 
** (updated 12/16/2010 -- program was not all good: 
**   * spurious "func" parameter in MakeGraph and "safetyarg" was missing.)
** ------------------------------------------------------------------------------------------ 
**  
** ============================================================================= */

/* =============================== Code: ==================================================== */
/* Complex function magnitude/phase plotter. */

/* To use:
*     1. Define function to graph as func(z).
*     2. Load this program.
*     3. Execute MakeGraph(width, height, x0, y0, x1, y1, filename) with the parameters given as follows:
*        width, height = width/height of image in pixels
*        x0, y0, x1, y1 = rectangle of complex plane to graph: x0 + y0i in upper-left corner to x1 + y1i in lower-right corner
*        filename = name of file to save as.
* Output is in .PPM format.
*/


/* Color conversion (HSB to RGB). */

HSB2RGB(hsb) = {
       local(H=hsb[1]);
       local(S=hsb[2]);
       local(B=hsb[3]);
       local(HH);
       local(F);
       local(P);
       local(Q);
       local(T);

       HH = floor(6*H)%6;
       F = (6*H) - floor(6*H);
       P = B*(1 - S);
       Q = B*(1 - (F*S));
       T = B*(1 - (1-F)*S);
       if(B > 1.0, B = 1.0);
       if(B < 0.0, B = 0.0);
       if(P > 1.0, P = 1.0);
       if(P < 0.0, P = 0.0);
       if(Q > 1.0, Q = 1.0);
       if(Q < 0.0, Q = 0.0);
       if(T > 1.0, T = 1.0);
       if(T < 0.0, T = 0.0);

       if(HH == 0, return([B, T, P]));
       if(HH == 1, return([Q, B, P]));
       if(HH == 2, return([P, B, T]));
       if(HH == 3, return([P, Q, B]));
       if(HH == 4, return([T, P, B]));
       if(HH == 5, return([B, P, Q]));
       }

/* Safe argument. */
safetyarg(z) = if(z == 0, 0, arg(z));

/* Make graph. */
MakeGraph(width, height, x0, y0, x1, y1, filename) = {
       xstep = (x1 - x0)/width;
       ystep = (y1 - y0)/height;
       write(filename, "P3");
       write(filename, "# ", filename);
       write(filename, width, " ", height);
       write(filename, "255");

       for(y=0, height-1,
           for(x=0, width-1,
                  xx = x0+(xstep*x);
                  yy = y0+(ystep*y);
               z = xx+yy*I;
               funcvalue = func(z);
               mag = abs(funcvalue);
               phase = safetyarg(funcvalue);
               H = phase/(2*Pi);
               S = 1/(1 + 0.3*log(mag + 1));
               B = 1 - 1/(1.1 + 5*log(mag + 1));
               RGB = HSB2RGB([H, S, B]);
                  Red = floor(RGB[1]*255.0);
                  Green = floor(RGB[2]*255.0);
                  Blue = floor(RGB[3]*255.0);
               write1(filename, Red, " ", Green, " ", Blue, "  ");
              );
           write(filename, "");
       );
       print("Done.");
    } 
