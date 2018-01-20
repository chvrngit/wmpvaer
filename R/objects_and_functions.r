#' p_kprs computes the first derivative as a function of s using an analytical formula 
#'
#' @param s The value needed for the first derivative to equal log(Wilks)
#' @param n The error Df of the one-way MANOVA analysis considered
#' @param m The hypothesis Df of the one-way MANOVA analysis considered
#' @param omega A vector of eigenvalues of the Wilks Non-Centrality Parameter corresponding
#'               to one independent variable.
p_kprs<-function(s,n,m,omega) {

a<-s;
b<-Rmpfr::mpfr( (((n+m)/2)+a), 120  );
x<-Rmpfr::roundMpfr(Rmpfr::mpfr(((-1/2)*omega),120  ),120  );
p<-length(x);

omega<-Rmpfr::mpfr(omega,120  );
q1a1<-Rmpfr::mpfr( ( ((omega-b)%*%(omega-b))%*%(matrix(1,1,p)) ),120  ) 
q1<-t( Rmpfr::mpfr(          as.matrix((4*(a)*omega)),           120  )  )
qm<-q1+q1a1;
sqm<-sqrt(qm);

dyds1a<-  (b - x) + sqm;
dyds1b<- (a)*((x+b)* ((sqm)^-1))   
dyds1c<-  dyds1a-a-dyds1b;
dyds1d<-  (dyds1a) ^ 2
dyds<-   2*(dyds1c/dyds1d)


yhat<-Rmpfr::mpfr((2*(a)*(dyds1a^(-1))), 120  )

cap_1<-((b -((p + 1)/4))*(p/b))+(p*log(b));

cap_m1<-Rmpfr::mpfr((matrix(0,p,p)),120  );
cap_m2<-Rmpfr::mpfr((matrix(0,p,p)),120  );
cap_m3<-Rmpfr::mpfr((matrix(0,p,p)),120  );
cap_m4<-Rmpfr::mpfr((matrix(0,p,p)),120  );

for(i in (1:p)){
for(j in (i:p)){
cap_m1[1,j]<-(  (  (yhat[1,i]*yhat[1,j])/a )  +   (    (1-yhat[1,i])*(1-yhat[1,j]) / (b-a) ) )^(-1);
}
}

for(i in (1:p)) {
for(j in (i:p)) {
cap_m2[i,j]<-(          (yhat[1,i]/a) -( (1-yhat[1,i])/(b-a) ) )*(dyds[1,j]);
}
}

for(i in (1:p)) {
for (j in (i:p)){
cap_m3[i,j]<- (  (yhat[1,j]/a)- ((1-yhat[1,j])/(b-a))  )*(dyds[1,i]);  
}
}

for(i in (1:p)) {
for (j in (i:p)) {
cap_m4[i,j]<-    (yhat[1,i]*yhat[1,j])/(a^2)   ;
}
}

cap_m5a<-(cap_m1)*(cap_m2+cap_m3-cap_m4);

#  then sum across columns

cap_m5b<-Rmpfr::mpfr((matrix(NA,p,1)),120  );
for(i in (1:p)) {
cap_m5b[i,1]<-Rmpfr::mpfr((sum(cap_m5a[i,(1:p)])),120  );
}

#  and sum across rows
cap_m5c<-Rmpfr::mpfr((sum(cap_m5b[(1:p),1])),120  );
# multiply by -1/2 and this is finished
cap_m5<-Rmpfr::mpfr((((-.5)*(cap_m5c))),120  );
#
cap_m61<-log(yhat/a);
cap_m62<-(1);
cap_m63a<-(a)*((yhat)^(-1));
cap_m63b<-(b-a)*((1-yhat)^(-1))
cap_m63<-(cap_m63a-cap_m63b)* dyds;
cap_m6a<-(cap_m61-cap_m62+cap_m63);


#initialize a place to store the answer 
# sum columns
cap_m6<-Rmpfr::mpfr((sum(cap_m6a[1,1:p])),120  );
# finally compute capdelta_s
capdelta_s<-Rmpfr::mpfr((cap_1+cap_m5+cap_m6),120  ); 
# digamma part
sc1<-Rmpfr::mpfr(((n/2)+a),120  );
sc2<-Rmpfr::mpfr((((n+m)/2)+a),120  );
ivec<-(1:p);
fac1<-Rmpfr::mpfr((digamma(sc1-((.5)*(ivec-1)) )),120  ); 
fac2<-Rmpfr::mpfr((digamma(sc2-((.5)*(ivec-1)) )),120  );
fac3<- (fac1-fac2);
digamma_part<-sum(fac3);
# put everything together
first_der<-digamma_part + capdelta_s;

return(first_der)
}


#' p_oneFhatone, computation of an approximate 1F1 value 
#'
#' @param s The value needed for the first derivative to equal log(Wilks) 
#' @param n The error Df of the one-way MANOVA analysis considered
#' @param m The hypothesis Df of the one-way MANOVA analysis considered
#' @param omega The a vector of eigenvalues of the Wilks Non-Centrality Parameter 
#'              corresponding to one independent variable.
p_oneFhatone<- function(s,n,m,omega){

a<-s;
b<-Rmpfr::mpfr(( ((n+m)/2)+a),precBits=120  );
x<-Rmpfr::roundMpfr(Rmpfr::mpfr(((-1/2)*omega),precBits=120  ),precBits=120  );

p<- length(x);

f1<-(b-x)+sqrt(((x-b)^2)+(4*a*x));
yht<-(2*a)*((f1)^(-1));
yhat<-t(as.matrix(yht));


m1<-matrix(.5,p,p);
M1 <-Rmpfr::mpfr((m1),precBits=120  )  ##-> "mpfrMatrix"


for(i in (1:p) ){
for(j in (i:p) ){
M1[i,j]<-(yhat[1,i] * (yhat[1,j]/a));
}
}


m2<-matrix(.5,p,p);
M2<-Rmpfr::mpfr((m2),precBits=120  )  ##-> "mpfrMatrix"


for(i in (1:p)) {
for(j in (i:p)) {
M2[i,j]<-( (1-yhat[1,i]) * ((1-yhat[1,j])/(b-a)) );
}
}

M3<-(M1+M2);

#second multipication

m4<-matrix(NA,p,1);
M4<-Rmpfr::mpfr((m4),precBits=120  )  ##-> "mpfrMatrix"

for (i in (1:p)){
M4[i,1]<-prod(M3[i,(1:p)]);
}

# final answer for R11

R11<-prod(M4[(1:p),1]);

#final computation for 1F1 

m5a<-(yhat/a)^(a);
m5b<-(1-yhat)/(b-a);
m5c<-exp(x*yhat);
m5<-(m5a*m5b*m5c) ;

pm5<-prod(m5);

oneF1a<-b^((p*b)-((p*(p+1)/4)));
oneF1b<-R11^(-(1/2));

oneF1<-oneF1a*oneF1b*pm5;

return(oneF1);

}

#' p_mgamma, computation of a multivariate gamma function needed by p_lmglm 
#'
#' @param p The size or length of a required eigenvector 
#' @param a A variable computed in p_lmglm
p_mgamma<- function(p,a)   {

ivec<-c(1:p);
ivecm<-Rmpfr::mpfr(ivec,120 );
ivecm1<-(1-ivecm)/2;
ivecm2<-a+ivecm1;

ivecm2<-t(as.matrix(ivecm2))

if(any(ivecm2<=0)==TRUE){
for(i in (1:ncol(ivecm2))){
if(ivecm2[1,i]<=0) {ivecm2[1,i]<-1} 
}
}

vg<-gamma(ivecm2)
mg<-prod(vg);
piv<-Rmpfr::mpfr(pi,120)
fct<-Rmpfr::mpfr((p*((p-1)/4)),120);
firstfact<-piv^fct;
g_u_pfa<- mg*firstfact;
g_u_pfa
return(g_u_pfa);
}

#' p_lmglm, computation of log(mgf) of the log(Wilks (Lambda) Statistic) 
#'
#' @param s_hat, The value needed for the first derivative to equal log(Wilks), found by function p_find_s
#' @param n The error Df of the one-way MANOVA analysis considered
#' @param m The hypothesis Df of the one-way MANOVA analysis considered
#' @param omega The a vector of eigenvalues of the Wilks Non-Centrality Parameter 
#'              corresponding to one independent variable.
p_lmglm<-function (s_hat,n,m,omega) {


if(is.na(s_hat)){
ans_lmglm<-NA;
} else {

p<-length(omega);
sv1<-Rmpfr::mpfr(((n/2)+s_hat),precBits=120  )
sv2<-Rmpfr::mpfr(((n+m)/2),precBits=120  )
sv3<-Rmpfr::mpfr((n/2),precBits=120  )
sv4<-Rmpfr::mpfr((sv2+s_hat),precBits=120  )
factor1<-p_mgamma(p,sv1);
factor2<-p_mgamma(p,sv2);
factor3<-p_mgamma(p,sv3);
factor4<-p_mgamma(p,sv4);
factor5<-p_oneFhatone(s_hat,n,m,omega);

ans_mglm1<-factor1*factor2;
ans_mglm2<-factor3*factor4;
ans_mglm3<-ans_mglm1/ans_mglm2;
ans_mglm<-abs(ans_mglm3*factor5);

ans_lmglm<-Rmpfr::mpfr((log(ans_mglm)),precBits=120  )

}
return(ans_lmglm);
}

#' find_min_s_hat finds minimum s_hat for bisection search,s such that First-Der(s)=log(Wilks)
#'
#' @param n Error DF for MANOVA
#' @param m The hypothesis Df of the one-way MANOVA analysis considered
#' @param omega The a vector of eigenvalues of the Wilks Non-Centrality Parameter corresponding 
#'              to one independent variable.
find_min_s_hat<-function(omega,n,m) {

s_hat_min<-Rmpfr::mpfr((-((n+m)/2)+.1),120  ); 
p<-length(omega);

b1<-Rmpfr::mpfr(((n+m)/2),120  )
b<-Rmpfr::mpfr((b1+s_hat_min),120  )
omega<-Rmpfr::mpfr(omega,120  );
q1a1<-Rmpfr::mpfr( ( ((omega-b)%*%(omega-b))%*%(matrix(1,1,p)) ),120  ) 
q1<-t( Rmpfr::mpfr(          as.matrix((4*(s_hat_min)*omega)),           120  )  )
qm<-q1+q1a1;


if (any(qm<=0)==TRUE) {
repeat {
s_hat_min<-Rmpfr::mpfr((s_hat_min+10),120  )
b<-Rmpfr::mpfr((b1+s_hat_min),120  );
q1a1<-Rmpfr::mpfr( ((omega-b)%*%(omega-b))%*%(matrix(1,1,p)),120  ) 
q1<-t(as.matrix(Rmpfr::mpfr((4*(s_hat_min)*omega),120  )))
qm<-q1+q1a1;
if (any(qm<=0)==FALSE) {break}
}
}

if (any(qm<=0)==FALSE) {
repeat{
s_hat_min<-Rmpfr::mpfr((s_hat_min-1),120  )
b<-Rmpfr::mpfr((b1+s_hat_min),120  );
q1a1<-Rmpfr::mpfr( ((omega-b)%*%(omega-b))%*%(matrix(1,1,p)),120  ) 
q1<-t(as.matrix(Rmpfr::mpfr((4*(s_hat_min)*omega),120  )))
qm<-q1+q1a1;
if (any(qm<=0)==TRUE) {break}
}
}

if (any(qm<=0)==TRUE) {
repeat{
s_hat_min<-Rmpfr::mpfr((s_hat_min+.1),120  ) #1
b<-Rmpfr::mpfr((b1+s_hat_min),120  );
q1a1<-Rmpfr::mpfr( ((omega-b)%*%(omega-b))%*%(matrix(1,1,p)),120  ) 
q1<-t(as.matrix(Rmpfr::mpfr((4*(s_hat_min)*omega),120  )))
qm<-q1+q1a1;
if (any(qm<=0)==FALSE) {break}
}
}

if (any(qm<=0)==FALSE) {
repeat{
s_hat_min<-Rmpfr::mpfr((s_hat_min-.01),120  ) #2
b<-Rmpfr::mpfr((b1+s_hat_min),120  );
q1a1<-Rmpfr::mpfr( ((omega-b)%*%(omega-b))%*%(matrix(1,1,p)),120  ) 
q1<-t(as.matrix(Rmpfr::mpfr((4*(s_hat_min)*omega),120  )))
qm<-q1+q1a1;
if (any(qm<=0)==TRUE) {break}
}
}

if (any(qm<=0)==TRUE) {
repeat{
s_hat_min<-Rmpfr::mpfr((s_hat_min+.001),120  ) #3
b<-Rmpfr::mpfr((b1+s_hat_min),120  );
q1a1<-Rmpfr::mpfr( ((omega-b)%*%(omega-b))%*%(matrix(1,1,p)),120  ) 
q1<-t(as.matrix(Rmpfr::mpfr((4*(s_hat_min)*omega),120  )))
qm<-q1+q1a1;
if (any(qm<=0)==FALSE) {break}
}
}


if (any(qm<=0)==FALSE) {
repeat{
s_hat_min<-Rmpfr::mpfr((s_hat_min-.0001),120  ) #4
b<-Rmpfr::mpfr((b1+s_hat_min),120  );
q1a1<-Rmpfr::mpfr( ((omega-b)%*%(omega-b))%*%(matrix(1,1,p)),120  ) 
q1<-t(as.matrix(Rmpfr::mpfr((4*(s_hat_min)*omega),120  )))
qm<-q1+q1a1;
if (any(qm<=0)==TRUE) {break}
}
}

if (any(qm<=0)==TRUE) {
repeat{
s_hat_min<-Rmpfr::mpfr((s_hat_min+.00001),120  ) #5
b<-Rmpfr::mpfr((b1+s_hat_min),120  );
q1a1<-Rmpfr::mpfr( ((omega-b)%*%(omega-b))%*%(matrix(1,1,p)),120  ) 
q1<-t(as.matrix(Rmpfr::mpfr((4*(s_hat_min)*omega),120  )))
qm<-q1+q1a1;
if (any(qm<=0)==FALSE) {break}
}
}


if (any(qm<=0)==FALSE) {
repeat{
s_hat_min<-Rmpfr::mpfr((s_hat_min-.000001),120  ) #6
b<-Rmpfr::mpfr((b1+s_hat_min),120  );
q1a1<-Rmpfr::mpfr( ((omega-b)%*%(omega-b))%*%(matrix(1,1,p)),120  ) 
q1<-t(as.matrix(Rmpfr::mpfr((4*(s_hat_min)*omega),120  )))
qm<-q1+q1a1;
if (any(qm<=0)==TRUE) {break}
}
}

if (any(qm<=0)==TRUE) {
repeat{
s_hat_min<-Rmpfr::mpfr((s_hat_min+.0000001),120  ) #7
b<-Rmpfr::mpfr((b1+s_hat_min),120  );
q1a1<-Rmpfr::mpfr( ((omega-b)%*%(omega-b))%*%(matrix(1,1,p)),120  ) 
q1<-t(as.matrix(Rmpfr::mpfr((4*(s_hat_min)*omega),120  )))
qm<-q1+q1a1;
if (any(qm<=0)==FALSE) {break}
}
}


if (any(qm<=0)==FALSE) {
repeat{
s_hat_min<-Rmpfr::mpfr((s_hat_min-.00000001),120  ) #8
b<-Rmpfr::mpfr((b1+s_hat_min),120  );
q1a1<-Rmpfr::mpfr( ((omega-b)%*%(omega-b))%*%(matrix(1,1,p)),120  ) 
q1<-t(as.matrix(Rmpfr::mpfr((4*(s_hat_min)*omega),120  )))
qm<-q1+q1a1;
if (any(qm<=0)==TRUE) {break}
}
}

if (any(qm<=0)==TRUE) {
repeat{
s_hat_min<-Rmpfr::mpfr((s_hat_min+.000000001),120  ) #9
b<-Rmpfr::mpfr((b1+s_hat_min),120  );
q1a1<-Rmpfr::mpfr( ((omega-b)%*%(omega-b))%*%(matrix(1,1,p)),120  ) 
q1<-t(as.matrix(Rmpfr::mpfr((4*(s_hat_min)*omega),120  )))
qm<-q1+q1a1;
if (any(qm<=0)==FALSE) {break}
}
}



if (any(qm<=0)==FALSE) {
repeat{
s_hat_min<-Rmpfr::mpfr((s_hat_min-.0000000001),120  ) #10
b<-Rmpfr::mpfr((b1+s_hat_min),120  );
q1a1<-Rmpfr::mpfr( ((omega-b)%*%(omega-b))%*%(matrix(1,1,p)),120  ) 
q1<-t(as.matrix(Rmpfr::mpfr((4*(s_hat_min)*omega),120  )))
qm<-q1+q1a1;
if (any(qm<=0)==TRUE) {break}
}
}

if (any(qm<=0)==TRUE) {
repeat{
s_hat_min<-Rmpfr::mpfr((s_hat_min+.00000000001),120  ) #11
b<-Rmpfr::mpfr((b1+s_hat_min),120  );
q1a1<-Rmpfr::mpfr( ((omega-b)%*%(omega-b))%*%(matrix(1,1,p)),120  ) 
q1<-t(as.matrix(Rmpfr::mpfr((4*(s_hat_min)*omega),120  )))
qm<-q1+q1a1;
if (any(qm<=0)==FALSE) {break}
}
}
return(s_hat_min)
}
#' p_find_s performs a bisection search, finds an s such that First-Der(s)=log(Wilks)
#'
#' @param wilks The value of the Wilks statistic for one independent variable
#' @param n The error Df of the one-way MANOVA analysis considered
#' @param m The hypothesis Df of the one-way MANOVA analysis considered
#' @param omega The a vector of eigenvalues of the Wilks Non-Centrality Parameter corresponding 
#'              to one independent variable.
p_find_s<-function(wilks,n,m,omega) {

wilksn<-Rmpfr::mpfr(wilks,120  )
p<-length(omega);

s_hat_min<-Rmpfr::mpfr((-((n+m)/2)+.1),120  );    

s_hat_max<- Rmpfr::mpfr(1000000,120);

b1<-Rmpfr::mpfr(((n+m)/2),120  )
b<-Rmpfr::mpfr((b1+s_hat_min),120  )
omega<-Rmpfr::mpfr(omega,120  );
q1a1<-Rmpfr::mpfr( ( ((omega-b)%*%(omega-b))%*%(matrix(1,1,p)) ),120  ) 
q1<-t( Rmpfr::mpfr(          as.matrix((4*(s_hat_min)*omega)),           120  )  )
qm<-q1+q1a1;

if (any(qm<=0)==TRUE) {
s_hat_min<-find_min_s_hat(omega,n,m) 
}

s_hat_mid<-Rmpfr::mpfr((((s_hat_min+s_hat_max)/2)+.1),120  );

# Bisection search section for analytical derivative using routine p_kprs

counter<-0;

difans1<-Rmpfr::mpfr((log(wilksn)),120  );

repeat{

counter<-(counter+1);

difans2<-p_kprs(s_hat_mid,n,m,omega);
difans<-difans1-difans2;

if ((difans<0)==TRUE ){
s_hat_max<- Rmpfr::mpfr(((s_hat_max+s_hat_mid)/2),120  );
s_hat_mid<-Rmpfr::mpfr((s_hat_min+ ((s_hat_max-s_hat_min)/2)),120  );
}

if ((difans>0)==TRUE ) {
s_hat_min<-Rmpfr::mpfr((s_hat_min+((s_hat_max-s_hat_mid)/2)),120  );
s_hat_mid<-Rmpfr::mpfr(((s_hat_max+s_hat_mid)/2),precBits=120  );
}


difmaxmin<-Rmpfr::mpfr((abs(s_hat_max-s_hat_min)),120  );
s_hat<-Rmpfr::mpfr(s_hat_mid,120  );

if (  ((abs(difans))==0)==TRUE){ break};
if (  ((abs(difans))<=.00000001)==TRUE){ break};
if (  (difmaxmin<=.000000000001)==TRUE ) {
s_hat<-Rmpfr::mpfr(NA,120  );
break}
if (  (counter==10000)==TRUE) {
s_hat<-Rmpfr::mpfr(NA,120  );
break}
}

ans<-Rmpfr::mpfr(s_hat,120  )

return(ans)
}

#' p_kprs2 computes the 2nd derivative as a function of s_hat, s found by the function p_find_s 
#'
#' @param s_hat The value needed for the first derivative to equal log(Wilks), s found by p_find_s
#' @param n The error Df of the one-way MANOVA analysis considered
#' @param m The hypothesis Df of the one-way MANOVA analysis considered
#' @param omega A vector of eigenvalues of the Wilks Non-Centrality Parameter 
#'              corresponding to one independent variable.

p_kprs2<- function(s_hat,m,n,omega) {


if(is.na(s_hat)){
second_der<-NA;
} else {


parma<-10^(-8);
parm1<-(s_hat+parma);
parm2<-(s_hat-parma);

first_der1<-(p_kprs(parm1,m,n,omega));
first_der2<-(p_kprs(parm2,m,n,omega));

diffr<-(first_der1)-(first_der2)

second_der<-  diffr/(2*parma)
 
return(second_der);

}

}

#' p_l_and_r computes the Lugannani and Rice approx. of power 
#'
#' @param wilks The value of the Wilks statistic for one independent variable
#' @param n The error Df of the one-way MANOVA analysis considered
#' @param m The hypothesis Df of the one-way MANOVA analysis considered
#' @param omega A vector of eigenvalues of the Wilks Non-Centrality Parameter 
#'              corresponding to one independent variable.
p_l_and_r<- function(wilks,n,m,omega) {

s_hat1<-p_find_s(wilks,n,m,omega);
s_hat2<-s_hat1[[1]];

s_hat<-Rmpfr::mpfr(s_hat2,precBits=120  );

if(is.na(s_hat)){
power_answers<-matrix(data=(c(NA,NA)),1,2);
} else {

wilksn<-Rmpfr::mpfr(wilks,precBits=120  );
y<-log(wilksn)
k_s_hat<-p_lmglm(s_hat,n,m,omega);
u1<-p_kprs2(s_hat,n,m,omega);
u<-(s_hat)*(sqrt(abs(u1)))
r1<-y*s_hat
r2<-r1-k_s_hat;
r3<-(2*(abs(r2)))^(1/2)
r<-(sign(s_hat))*r3;
mpower1<-Rmpfr::mpfr( (Rmpfr::pnorm(r,lower.tail=FALSE)) ,precBits=120    );
piv<-Rmpfr::mpfr(pi,120)
nd_r1<-Rmpfr::mpfr((1/(sqrt(2*piv))),precBits=120  )
nd_r2<-(exp(-(r^2)/2) );
mpower2<-nd_r1*nd_r2
mpower3<-( ((1/u)-(1/r)))  ; 
mpower<-mpower1+(mpower2*mpower3); 
type2error<-(1-(mpower));
if((mpower>1)==TRUE){
mpower<-1;
type2error<-0;
}
if((type2error>1)==TRUE){
mpower<-0;
type2error<-1;
}
ans<-c(mpower,type2error);
power_answers<-t(as.matrix(ans));
}

return(power_answers)
}

#' Wilks MANOVA p-Value Alternative Estimation  
#'
#' \code{wmpvaerf} takes Wilks values from R's MANOVA function, adds power estimates and alternative p-values 
#' @param formula Same formula as used with the R MANOVA function
#' @section Warning:
#' This code example,on a single core, can take 30 seconds or more to complete. The program exploits two or more cores
#' if available,however.On a three core Intel system the code example takes about 11 seconds.  
#' @examples
#'     suppressPackageStartupMessages(library(Rmpfr)); 
#'     library(doParallel);
#'     library(wmpvaer);
#'     data(plastic);
#'     Y<-as.matrix(plastic[,1:3]);
#'     rate<-matrix(data=(plastic$rate),nrow(plastic),1);
#'     additive<-matrix(data=(plastic$additive),nrow(plastic),1); 
#'     wmpvaerf(Y~ rate*additive);
#' @keywords multivariate    &   Multivariate Techniques
#' @return list with three elements (matricies) 
#' @export
wmpvaerf<-function(formula) {
fit <- stats::manova(formula); 
eigs<-summary(fit)$Eigenvalues;
wstat<-summary(fit, test="Wilks")$stats;
leigs<-nrow(as.matrix(eigs));
weigs<-ncol(as.matrix(eigs));
dfm<-as.matrix(wstat[,1]);
wilksstat<-as.matrix(wstat[,2]);
ldf<-nrow(dfm);
n<-dfm[ldf,];
m<-(as.matrix(dfm[1:leigs,]));
w<-Rmpfr::asNumeric((as.matrix(wilksstat[1:leigs,])));
pmatrix<-as.matrix(wstat[1:leigs,6]);
lnw<-Rmpfr::asNumeric(log(w));
wilksv<-cbind(w,lnw)
rvec<-matrix(data=c("         ","         "),1,2);
wilksva<-Rmpfr::rbind(wilksv,rvec);
rownames(wilksva)<-rownames(wstat);
wilksvb<-Rmpfr::cbind(dfm,wilksva)
colnames(wilksvb)<-c("Df"," Wilks"," ln(Wilks)");
wilksvc<-noquote(wilksvb)

#section using parallel computation

n_usable_cores<-(parallel::detectCores( ))-1;


'%dopar%'<-'foreach::%dopar%'

if(n_usable_cores<2){

#create dummy matrix to hold results
resultsa<-Rmpfr::mpfr(matrix(NA,leigs,2),precBits= 120  );

#fill matrix
for (i in 1:leigs){
ans<-p_l_and_r(w[i,1],n,m[i,1],eigs[i,]);
resultsa[i,1:2]<-Rmpfr::mpfr(ans,precBits= 120  );
}

} else {

mincl<-min(leigs,n_usable_cores);
cl<-parallel::makeCluster(mincl);
doParallel::registerDoParallel(cl);
resultsa<-{
foreach::foreach(i=1:leigs,.combine=Rmpfr::rbind,
.export=c("p_find_s","p_kprs", "p_kprs2","p_mgamma","p_lmglm",
"p_l_and_r","p_oneFhatone" ) )  %dopar% 
{(p_l_and_r(w[i,1],n,m[i,1],eigs[i,]));}
}

parallel::stopCluster(cl);

}

# first matrix

if (is.numeric(leigs)==TRUE) { 
cvec<-matrix("CHR",leigs,1);}

resultsb<-Rmpfr::cbind(resultsa,cvec);
resultsc<-matrix(data=resultsb[,1:2],nrow(resultsa),2)
resultsd<-as.vector(resultsc);
resultsd<-replace(resultsd,resultsd=="NaN","No Answer")

resultse<-matrix(data=resultsd, nrow=leigs,ncol=2);
rownames(resultse)<-rownames(eigs);
colnames(resultse)<-c("Pr(x>ln(Wilks)),Power,precBits=120","Pr(x<=ln(Wilks)),p Value,precBits=120");

resultsf1<-format(resultse,justify="right")
resultsf<-noquote(resultsf1)

# second matrix
resultsa1a<-Rmpfr::cbind(resultsa,pmatrix)
resultsa1<-Rmpfr::format(resultsa1a, digits=7)
resultsd1<-as.vector(resultsa1);
resultsd1<-replace(resultsd1,resultsd=="NaN","No Answer")

resultse1<-matrix(data=resultsd1, nrow=leigs,ncol=3);
rownames(resultse1)<-rownames(eigs);
colnames(resultse1)<-c("Power,formatted"," p,formatted","  r-manova~p");
resultsf1a<-format(resultse1,justify="right")
resultsf1b<-noquote(resultsf1a)

answer<-list(wilksvc,resultsf,resultsf1b)
return(answer);
}

#' computeWilks takes a matrix of eignvalues and computes values of the Wilks statistic 
#'
#' @param omega Eigenvalues of an (E-inv)*H matrix,MANOVA analysis, Wilks non-centrality
#'   parameters
computeWilks<-function(omega){

p<-nrow(omega);
cols<-ncol(omega);

omega<-Rmpfr::roundMpfr(Rmpfr::mpfr(omega,120),120 );
omegap<-1/(1+omega);

omegam1<-matrix(NA,p,1);
wilksvs<-Rmpfr::mpfr((omegam1),precBits=120  )  ##-> "mpfrMatrix"

for (i in (1:p)){
wilksvs[i,1]<-prod(omegap[i,(1:cols)]);
}
return(wilksvs);
}

#' valid_mlandr contains validity tests for the wmpvaer::mlandr object 
#'
#' @param object The object in question is the S4 wmpvaer::mlandr
#'  
valid_mlandr<- function(object) {

m_df<-as.numeric(object@m_df);
eigs<-as.matrix(object@eigs);
eigs_num<-as.numeric(eigs);
e_df<-as.numeric(object@e_df);

msg <- NULL;
valid <- TRUE;

if( (any(is.na(eigs_num))) ==TRUE ){
 valid <- FALSE
 msg <- c(msg,"no characters allowed in in eigenvalue matrix") 
}

if (length(m_df) != nrow(eigs)) {
 valid <- FALSE
 msg <- c(msg,"Length model DF vector must equal number rows in eigs matrix");
}
 
if (any(m_df>e_df)==TRUE) {
 valid <- FALSE
 msg <- c(msg,"Model DF values must be < Error DF");
}

if ( (ncol(eigs)<2)==TRUE) {
 valid <- FALSE
 msg <- c(msg,"Analysis must have 2 or more dependent variables");
}

if (valid==FALSE) {print(msg)};
}

#' mlandr  S4 object 
#'
#' \code{mlandr} the S4 object takes a matrix of eignvalues, a vector of model DF values and a numeric value of error DF
#'        and returns Wilks statistics & BW power and p values 
#' @param e_df error DF for a MANOVA model, numeric vector , or 1x1 numeric matrix
#' @param m_df model DFs for a MANOVA model, a numeric vector,length >=1, or a matrix ncol=1 & nrow>=1,
#' @param eigs Eigenvalues of an (E-inv)*H matrix,MANOVA analysis, Wilks non-centrality parameter, matrix w. ncol>=2
mlandr<-methods::setClass("mlandr", slots=list(e_df="numeric",m_df="numeric",eigs="matrix") );
methods::setValidity("mlandr", valid_mlandr)
methods::setGeneric(name="compute.pwilks",def=function(object){standardGeneric("compute.pwilks")})
methods::setMethod(f="compute.pwilks","mlandr",
definition= function(object) {

eigs<-as.matrix(object@eigs);
w<-computeWilks(eigs);
leigs<-nrow(as.matrix(object@eigs));
weigs<-ncol(as.matrix(object@eigs));
m_df<-as.matrix(object@m_df);
e_df<-as.matrix(object@e_df);
dfm<-rbind(m_df,e_df);
ldf<-nrow(dfm);
n<-dfm[ldf,];
m<-(as.matrix(dfm[1:leigs,]));

lnw<-Rmpfr::mpfr((log(w)),120);
wilksv<-Rmpfr::cbind(w,lnw);
wilksv<-Rmpfr::format(wilksv,digits=9);
rvec<-matrix(data=c("         ","         "),1,2);
wilksva1<-Rmpfr::rbind(wilksv,rvec);
wilksva<-Rmpfr::cbind(dfm,wilksva1);
wilksvc<-noquote(wilksva);
colnames(wilksvc)<-c("Df","Wilks","ln(Wilks)");
rownames(wilksvc)<-c(rownames(eigs),"Residuals")

#section using parallel computation


n_usable_cores<-(parallel::detectCores( ))-1;

'%dopar%'<-'foreach::%dopar%'

if(n_usable_cores<2){

#create dummy matrix to hold results
resultsa<-Rmpfr::mpfr(matrix(NA,leigs,2),precBits= 120  );

#fill matrix
for (i in 1:leigs){
ans<-p_l_and_r(w[i,1],n,m[i,1],eigs[i,]);
resultsa[i,1:2]<-Rmpfr::mpfr(ans,precBits= 120  );
}

} else {

mincl<-min(leigs,n_usable_cores);
cl<-parallel::makeCluster(mincl);
doParallel::registerDoParallel(cl);
resultsa<-{
foreach::foreach(i=1:leigs,.combine=Rmpfr::rbind,
.export=c("p_find_s","p_kprs", "p_kprs2","p_mgamma","p_lmglm",
"p_l_and_r","p_oneFhatone" ) )  %dopar% {(p_l_and_r(w[i,1],n,m[i,1],eigs[i,]));}
}

parallel::stopCluster(cl);

}

# first matrix

if (is.numeric(leigs)==TRUE) { 
cvec<-matrix("CHR",leigs,1);}

resultsb<-Rmpfr::cbind(resultsa,cvec);
resultsc<-matrix(data=resultsb[,1:2],nrow(resultsa),2)
resultsd<-as.vector(resultsc);
resultsd<-replace(resultsd,resultsd=="NaN","No Answer");

resultse<-matrix(data=resultsd, nrow=leigs,ncol=2);
rownames(resultse)<-rownames(eigs);
colnames(resultse)<-c("     Pr(x>ln(Wilks)),Power,precBits=120","     Pr(x<=ln(Wilks)),p Value,precBits=120");

resultsf1<-format(resultse,justify="right")
resultsf<-noquote(resultsf1)

# second matrix
resultsa1<-Rmpfr::format(resultsa, digits=7)
resultsd1<-as.vector(resultsa1);
resultsd1<-replace(resultsd1,resultsd=="NaN","No Answer")

resultse1<-matrix(data=resultsd1, nrow=leigs,ncol=2);
rownames(resultse1)<-rownames(eigs);
colnames(resultse1)<-c("Power,formatted"," p,formatted");
resultsf1a<-format(resultse1,justify="right")
resultsf1b<-noquote(resultsf1a)

resultsf1c<-Rmpfr::rbind(resultsf1a,rvec);


wilksvd<-Rmpfr::cbind(wilksvc,resultsf1c[,2]);
colnames(wilksvd)<-c("Df","   Wilks","   ln(Wilks)"," p,formatted");
rownames(wilksvd)<-c(rownames(eigs),"Residuals");
wilksvd<-noquote(wilksvd)

answer<-list(wilksvd,resultsf)
return(answer);
});


#' mlandrf  mlandr S4 object R function wrapper 
#'
#' \code{mlandrf} takes a matrix of eignvalues, a vector of model DF values and a numeric value of error DF
#'        and returns Wilks statistics & BW power and p values by wrapping the S4 object 
#'        wmpvaer::mlandr in an R function.
#' @param e_df error DF for a MANOVA model, numeric vector , or 1x1 numeric matrix
#' @param m_df model DFs for a MANOVA model, a numeric vector,length >=1, or a matrix ncol=1 & nrow>=1,
#' @param eigs Eigenvalues of an (E-inv)*H matrix,MANOVA analysis, Wilks non-centrality parameter, matrix w. ncol>=2
#' @examples
#'     suppressPackageStartupMessages(library(Rmpfr)); 
#'     library(doParallel);
#'     library(wmpvaer);
#'     library(methods);
#'     data(eigs_plastic);
#'     eigsm<-as.matrix(eigs_plastic)
#'     colnames(eigsm)<-c(" "," "," ");
#'     rownames(eigsm)<-c("rate","additive","rate:additive");
#'     e_dfm<-16;
#'     m_dfm<-c(1,1,1);
#'     mlandrf(e_dfm,m_dfm,eigsm);
#' @keywords multivariate    &   Multivariate Techniques
#' @return list with two elements (matricies) 
#' @export
mlandrf<-function(e_dfm,m_dfm,eigsm) {
answer<-compute.pwilks(new("mlandr",e_df=e_dfm,m_df=m_dfm,eigs=eigsm));
return(answer);
}

#' oneway_mss   Computes est (based on the Wilks statistic) of # replications and sample 
#'               estimates for 1-way MANOVA w. 2 to 7 dependent vars and from 2 to 12 
#'               levels of an independent var.
#' \code{oneway_mss} takes two parameters: a scalar indicating number of levels and 
#'                 a vector of the eigenvalues of a 1-way MANOVA (E-Inv)*H matrix. 
#' @param n_level Is a numeric scalar, =>2 & =<12 
#' @param omega Is numeric vector,eigenvalues of a 1-way MANOVA (E-Inv)*H matrix, length=>2 & =<7 
#' @examples
#'     suppressPackageStartupMessages(library(Rmpfr)); 
#'     library(wmpvaer);
#'     rate_eigs<-c(1.6187719, -9.273797e-17,  3.365044e-18);
#'     rate_n_level<-2;  
#'     oneway_mss(rate_n_level,rate_eigs);
#' @keywords multivariate    &   Multivariate Techniques
#' @return (if replications<=20) printed results. 
#' @export
oneway_mss<- function(n_level,omega) {

if((is.atomic(n_level) && (length(n_level) == 1L))==FALSE){stop("n_level must be a scalar")};

if( (is.numeric(n_level))==FALSE){stop("n_level must be numeric") };

omega<-as.numeric(omega);
if( (is.vector(omega))==FALSE ){stop("omega must be a vector")};
if( (any(is.na(omega)))==TRUE ){stop("omega must be a numeric")};


if( (n_level<2)==TRUE ){stop("1-way levels must >=2") };
if( (n_level>12)==TRUE ){stop("program limited to levels <=12") };

l_omega<-length(omega);

if( (l_omega<2)==TRUE ){stop("model must have 2 or more dependent vars") };
if( (l_omega>7)==TRUE ){stop("program restricted to models between 2 and 7 dependent vars") };

m<-(n_level-1);

w_mid<-wmpvaer::alpham[m,(l_omega+2)]

J<-1;
repeat{
J<-J+1;
if((J>20)==TRUE){
cat("\n")
cat("1-Way MANOVA Balanced A-priori Study Design","\n")
cat("----------------------------","\n")
cat("eigs(E-inv*H)=",omega,"\n")
cat("levels=",n_level,",dependent vars=",l_omega,"\n")
cat("alpha approx.=.05",",power>=.9","\n")
cat("replications needed=NO ANSWER","\n")
cat("Program Stopped, est. reps>20")
break}

ndf<-(n_level*J)

difans2<-p_l_and_r(w_mid,ndf,m,((J*omega)))[1,1];

if (  (difans2>=.9)==TRUE ){
cat("\n")
cat("1-Way MANOVA Balanced A-priori Study Design","\n")
cat("----------------------------","\n")
cat("eigs(E-inv*H)=",omega,"\n")
cat("levels=",n_level,",dependent vars=",l_omega,"\n")
cat("alpha approx.=.05",",power>=.9","\n")
cat("replications needed=",J,"\n")
cat("total sample size=",(n_level*J));
#answer<-c( J, ndf)
#return(answer)
break}
}


}

