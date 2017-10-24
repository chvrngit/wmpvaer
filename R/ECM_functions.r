#' p_kprsECM computes the first derivative as a function of s using an analytical formula 
#' @param s The value needed for the first derivative to equal log(lambdaECM)
#' @param N1 Number of rows in matrix 1  
#' @param N2 Number of rows in matrix 2  
#' @param sigma_eigv A vector of eigenvalues of (Sigma11)*(inv-Sigma22)
p_kprsECM<-function(s,N1,N2,sigma_eigv) {

p<-length(sigma_eigv);
n1<-N1-1;
n2<-N2-1;
n<-n1+n2;

s<-Rmpfr::mpfr(s,120)
sigma_eigv<-Rmpfr::mpfr(sigma_eigv,precBits=120  );
onemsig<-(1-sigma_eigv);

x<-Rmpfr::roundMpfr(onemsig,120);
x<-t(as.matrix(x));

#computation of the derivative of the log first term

 .f1 <- 2 * (s/n); 
 .f4 <- n * (1 + .f1)/2; 
 .f6 <- .f4 - (1 + p)/4; 
 .f7 <- p * .f6;
 .f8 <- .f4^.f7; 
 
cap_1<- (p*((log(n) + log1p(.f1) - 0.693147180559945)*.f8 + .f4^(.f7- 1)*.f6)/.f8) ;


#computation of the derivative of log(R21)

cap_m1<-Rmpfr::mpfr((matrix(0,p,p)),120  );

for(i in (1:p)){
for(j in (i:p)){

 .a1 <- 1 + 2 * (s/n);
 .a2 <- n/2;
 .a5 <- n1 * .a1/2 - s; 
 .a6 <- x[1,i];
 .a7 <- x[1,j]; 
 .a9 <- n * .a1/2;
 .a10 <- .a2 - n1/2; 
 .a11 <- .a5 * .a7; 
 .a12 <- .a5 * .a6;
 .a13 <- s * .a10; 
 .a14 <- .a1*(.a2-4*(.a13*.a6));
 .a15 <- .a11 - .a9; 
 .a16 <- .a12 - .a9;
 .a17 <- n1/n; 
 .a21 <- .a15^2; 
 .a22 <- sqrt(abs(.a16^2 + .a14 - .a12));
 .a25 <- sqrt(abs(.a21 + .a1 * (.a2 - 4 * (.a13 * .a7)) - .a11));
 .a26 <- .a17 - 1;
 .a32 <- 4 * (.a1 * .a10) + 4 * (s * (1 - .a17));
 .a33 <- sqrt(abs(.a21 + .a14 - .a11));
 .a34 <- .a26 * .a7; 
 .a35 <- 1 - .a22;
 .a36 <- .a9 - s;
 .a38 <- .a32 + .a17 - 1; 
 .a40 <- 1 - .a25;
 .a41 <- 1 + 2 * (.a15 * (.a34 - 1)); 
 .a44 <- 1 - .a22 * .a6;
 .a45 <- 1 - .a25 * .a7;
 .a47 <- 1 + 2 * (.a16 * (.a26 * .a6 - 1)) - .a38 * .a6; 
 .a51 <- .a41 - (.a32 * .a6 + .a34);
 .a52 <- .a41 - .a38 * .a7;
 .a53 <- s * .a36;

cap_m1[i,j]<-(((0.5*(.a51/.a33)-.a33/s)*.a22+0.5*(.a47*.a33/.a22))/s-((0.5*(.a35*.a51/.a33)+0.5*((1-.a33)*.a47/.a22))/.a36+n1*
 (((((0.5*(.a6/.a44)+0.5/.a22)*.a35-0.5)*.a47*.a25+0.5*(.a35 * .a52 * .a22/.a25))*.a40+(0.5*(.a40*.a7/.a45)-0.5)*.a35*.a52*.a22)*
 .a1/(2*.a53)+.a35*.a40*(1/(n*s*.a36)-.a1*.a36/(2*.a53^2))*.a22*.a25)*.a6*.a7/(.a44*.a45)))/
 ((1-(n1*.a40*.a1*.a22*.a25*.a6*.a7/(2*(s*.a44*.a45))+.a33))*.a35/.a36+.a22*.a33/s  ) ;

}
}

#  then sum across columns

cap_m5b<-Rmpfr::mpfr((matrix(NA,p,1)),120  );
for(i in (1:p)) {
cap_m5b[i,1]<-Rmpfr::mpfr((sum(cap_m1[i,(1:p)])),120  );
}

#  and sum across rows
cap_m5c<-Rmpfr::mpfr((sum(cap_m5b[(1:p),1])),120  );
# multiply by -1/2 and this is finished
cap_m5<-Rmpfr::mpfr((((-.5)*(cap_m5c))),120  );
#

thrdm<-Rmpfr::mpfr((matrix(NA,1,p)),120 );

for(i in (1:p)){

.c1 <- 1 + 2 * (s/n);
.c2 <- x[1,i]; 
.c4 <- (n1 * .c1/2 - s) * .c2 - n * .c1/2;
.c6 <- .c4^2/s; 
.c7 <- .c6^(s/2); 

thrd1<-((0.5*(.c6^((s - 1)/2)*.c4*(2*((n1/n - 1)*.c2 - 1)-.c4/s)/sqrt(.c6)) + 0.5 * (.c7 * (2 * log(abs(.c4/s)))))/.c7);
#thrd1<-((0.5*(.c6^((s - 1)/2)*.c4*(2*((n1/n - 1)*.c2 - 1)-.c4/s)/sqrt(.c6)))/.c7);

.d1 <- 1 + 2 * (s/n);
.d2 <- x[1,i]; 
.d3 <- (n1 * .d1/2 - s) * .d2; 
.d4 <- n/2; 
.d6 <- n * .d1/2;
.d7 <- .d3 - .d6;
.d8 <- .d4 - n1/2; 
.d9 <- n1/n; 
.d13 <- sqrt(abs((.d7^2) +.d1*(.d4-4*(s*.d8*.d2)) - .d3));

thrd2<- (-(0.5*((1 + 2*(.d7 *((.d9 - 1)*.d2-1))-(4 * (.d1*.d8) + 4 * (s*(1-.d9))+.d9- 1)*.d2)*(.d6-s)/((1-.d13)*.d13)))); 

.e1 <- 1 + 2 * (s/n);
.e2 <- x[1,i];
.e4 <- n1 * .e1/2; 
.e5 <- (.e4 - s) * .e2;
.e6 <- n/2; 
.e7 <- .e5 - n * .e1/2;
.e8 <- .e6 - n1/2;
.e12 <- sqrt(abs(.e7^2 + .e1 * (.e6 - 4 * (s * .e8 * .e2)) - .e5));
.e13 <- abs(1 - .e12 * .e2); 
.e14 <- n1/n;
.e15 <- .e13^.e4;

thrd3<-(n1*((1+2*(.e7*((.e14-1)*.e2-1))-(4*(.e1*.e8)+4*(s*(1-.e14))+.e14-1)*.e2)*.e1*.e2/(4*(.e13^(1+.e4)*.e12))-log(.e13)/(n*.e15))*.e15);

thrdm[1,i]<-thrd2+thrd3;
}

# sum columns
cap_m6<-Rmpfr::mpfr((sum(thrdm[1,1:p])),120  );
# finally compute capdelta_s
capdelta_s<-Rmpfr::mpfr((cap_1+cap_m5+cap_m6),120  ); 
#capdelta_s<-Rmpfr::mpfr((cap_1+cap_m6),120  ); 


# digamma part
dgm<-Rmpfr::mpfr((matrix(NA,1,p)),120 );
for(i in (1:p)){
.g1 <- (i - 1)/2;
.g2 <- 1 + 2 * (s/n); 
dgm[1,i]<-(n1*digamma(n1*(.g2/2)-.g1) + n2*digamma(n2*(.g2/2)-.g1))/(n- digamma(n*(.g2/2)-.g1))  
}
digamma_part<-sum(dgm[1,1:p])

deltadet<-prod(sigma_eigv);
deltadet_part<-n1*log(deltadet)/n;

# put everything together
first_der<-digamma_part + capdelta_s+deltadet_part;

return(first_der)
}
#' p_twoF1ECM, computation of an approximate 2F1 value for equality of covariances 
#' @param s The value needed for the first derivative to equal log(lambdaECM) 
#' @param N1 number of rows in data matrix N1
#' @param N2 number of rows in data matrix N2
#' @param sigma_eigv A vector of eigenvalues of (Sigma11)*(inv-Sigma22)
p_twoF1ECM<-function(s,N1,N2,sigma_eigv) {

p<-length(sigma_eigv);
n1<-N1-1;
n2<-N2-1;
n<-n1+n2;

s<-Rmpfr::mpfr(s,120)

sigma_eigv<-Rmpfr::mpfr(sigma_eigv,precBits=120  );
onemsig<-abs(1-sigma_eigv);

a<-s;
b<-((n1/2)*(1+((2*s)/n)));
c<-((n/2)*(1+((2*s)/n)));
mfact<-((b)/(a*(c-a)));

x<-Rmpfr::roundMpfr(onemsig,120);
x<-t(as.matrix(x))


#yhat calculations

tau<-(x*(b-c))-c;
f1<-sqrt((tau^2)- ((4*a*x)*(c-b)))- tau;
yht<-(2*a)*((f1)^(-1));
yhat<-(as.matrix(yht));


#R21 calculations

m1<-matrix((1/3),p,p);
M1<-Rmpfr::mpfr((m1),precBits=120  )  ##-> "mpfrMatrix"

for(i in (1:p) ){
for(j in (i:p) ){
M1[i,j]<-(yhat[1,i] * (yhat[1,j]/a));
}
}

m2<-matrix((1/3),p,p);
M2<-Rmpfr::mpfr((m2),precBits=120  )  ##-> "mpfrMatrix"


for(i in (1:p)) {
for(j in (i:p)) {
M2[i,j]<-( (1-yhat[1,i]) * ((1-yhat[1,j])/(c-a)) );
}
}

m3<-matrix((1/3),p,p);
M3<-Rmpfr::mpfr((m3),precBits=120  )  ##-> "mpfrMatrix"


for(i in (1:p)) {
for(j in (i:p)) {


M3[i,j]<-((-1)*mfact*( (x[1,i]*yhat[1,i]*(1-yhat[1,j])) /
(1-(x[1,i])*yhat[1,j]))*( ( (x[1,j]*yhat[1,j]*(1-yhat[1,i])) /
(1-(x[1,j])*yhat[1,i]))) );

}
}


M4<-(M1+M2+M3);

# multipication

m5<-matrix(NA,p,1);
M5<-Rmpfr::mpfr((m5),precBits=120  )  ##-> "mpfrMatrix"

for (i in (1:p)){
M5[i,1]<-abs(prod(M4[i,(1:p)]));
}

# final answer for R21
R21<-(prod(M5[(1:p),1]));

#final computation for 2F1 
m6a<-(yhat/a)^(a);
m6b<-((1-yhat)/(c-a))^(c-a);
m6c<-((1-(x*yhat))^(-b));
m6<-(m6a*m6b*m6c);
pm6<-prod(m6);

twoF1a<-c^((p*c)-((p*(p+1)/4)));
twoF1b<-R21^(-(1/2));
twoF1c<-twoF1a*twoF1b*pm6;

#additional adjustment required
deltadet<-(prod(sigma_eigv))^((n1*s)/n)

twoF1<-deltadet*twoF1c;
return(twoF1);
}

#' p_mgammaECM, computation of a multivariate gamma function needed by p_lmglmECM 
#' @param p The size or length of a required eigenvector 
#' @param a A variable computed in p_lmglmECM
p_mgammaECM<- function(p,a)   {

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

#' p_lmglmECM, computation of log(mgf) of the log(lambdaECM) statistic 
#' @param s_hat, The value needed for the first derivative to equal log(lambdaECM), found by function p_find_sECM
#' @param N1 Number of rows in matrix 1  
#' @param N2 Number of rows in matrix 2  
#' @param sigma_eigv A vector of eigenvalues of (Sigma11)*(inv-Sigma22)
p_lmglmECM<-function(s_hat,N1,N2,sigma_eigv){

if(is.na(s_hat)){
ans_lmglm<-NA;
} else {

p<-length(sigma_eigv);
n1<-N1-1;
n2<-N2-1;
n<-n1+n2;

s_hat<-Rmpfr::mpfr(s_hat,120)
sigma_eigv<-Rmpfr::mpfr(sigma_eigv,precBits=120  );

mf<-(1+((2*s_hat)/n))

sv1<-Rmpfr::mpfr((n/2),precBits=120 );
sv2<-Rmpfr::mpfr(((n1/2)*mf),precBits=120 );
sv3<-Rmpfr::mpfr(((n2/2)*mf),precBits=120 );
sv4<-Rmpfr::mpfr((n1/2),precBits=120 );
sv5<-Rmpfr::mpfr((n2/2),precBits=120 );
sv6<-Rmpfr::mpfr(((n/2)*mf),precBits=120 );

factor1<-p_mgammaECM(p,sv1);
factor2<-p_mgammaECM(p,sv2);
factor3<-p_mgammaECM(p,sv3);
factor4<-p_mgammaECM(p,sv4);
factor5<-p_mgammaECM(p,sv5);
factor6<-p_mgammaECM(p,sv6);


factor7<-p_twoF1ECM(s_hat,N1,N2,sigma_eigv);

ans_mglm1<-factor1*factor2*factor3;
ans_mglm2<-factor4*factor5*factor6;
ans_mglm3<-ans_mglm1/ans_mglm2;
ans_mglm<-abs(ans_mglm3*factor7);

ans_lmglm<-Rmpfr::mpfr((log(ans_mglm)),precBits=120  )

}
return(ans_lmglm);
}
#' p_find_sECM performs a bisection search, finds an s such that FirstDer(s)=log(lambdaECM)
#' @param lambdaECM The value of the ECM statistic. 
#' @param N1 Number of rows in matrix 1  
#' @param N2 Number of rows in matrix 2  
#' @param sigma_eigv A vector of eigenvalues of (Sigma11)*(inv-Sigma22)
p_find_sECM<-function(lambdaECM,N1,N2,sigma_eigv) {

p<-length(sigma_eigv);
n1<-N1-1;
n2<-N2-1;
n<-n1+n2;

sigma_eigv<-Rmpfr::mpfr(sigma_eigv,precBits=120  );

lambdaECMn<-Rmpfr::mpfr(lambdaECM,120)

difans1<-Rmpfr::mpfr((log(lambdaECMn)),120  );
difansmx1<-Rmpfr::mpfr(matrix(NA,1,2),120)

counter_outer<-0;
difans1<-Rmpfr::mpfr((log(lambdaECM)),120  );

repeat{

set.seed( abs(as.numeric(Sys.time( ))-as.numeric(Sys.Date( ))));
searchsample<-runif(200,min=-10,max=10);

counter_outer<-counter_outer+1;

for(s_hat in searchsample){ 

difans2<-Rmpfr::mpfr(p_kprsECM(s_hat,N1,N2,sigma_eigv),120);

if((is.nan(difans2))==TRUE){next}
if((is.infinite(difans2))==TRUE){next}
difans<-Rmpfr::mpfr((difans1-difans2),120);
difansmx2<-Rmpfr::mpfr(matrix(NA,1,2),120)
difansmx2[1,2]<-abs(difans);
difansmx2[1,1]<-s_hat;
difansmx1<-Rmpfr::rbind(difansmx1,difansmx2)
}

difansmx1<-difansmx1[-1,]
md<-min(abs(difansmx1[,2]))

reduce<-difansmx1[abs(difansmx1[,2])==md,,drop=FALSE]
s_hat_mid<-reduce[1,1];

smm<-c((.95*s_hat_mid),(1.05*s_hat_mid))
s_hat_max<-max(smm);
s_hat_min<-min(smm);

difmaxmin_m<-Rmpfr::mpfr((matrix(NA,10000,1)),120 );

counter_inner<-0
difans_m<-Rmpfr::mpfr(matrix(NA,1000,1),120)


repeat{

counter_inner<-(counter_inner+1);

difans2<-p_kprsECM(s_hat_mid,N1,N2,sigma_eigv);
difans<-difans1-difans2;
difans_m[counter_inner,1]<-difans

if(counter_inner>1){
if ((abs(difans_m[counter_inner,1]-difans_m[counter_inner-1,1])<.0000000001)){
s_hat<-Rmpfr::mpfr(NA,120  );
break}
}

if ((difans<0)==TRUE ){
s_hat_max<- Rmpfr::mpfr(((s_hat_max+s_hat_mid)/2),120  );
s_hat_mid<-Rmpfr::mpfr((s_hat_min+ ((s_hat_max-s_hat_min)/2)),120  );
}

if ((difans>0)==TRUE ) {
s_hat_min<-Rmpfr::mpfr((s_hat_min+((s_hat_max-s_hat_mid)/2)),120  );
s_hat_mid<-Rmpfr::mpfr(((s_hat_max+s_hat_mid)/2),precBits=120  );
}


difmaxmin<-Rmpfr::mpfr((abs(s_hat_max-s_hat_min)),120  );
difmaxmin_m[counter_inner,1]<-difmaxmin;

s_hat<-Rmpfr::mpfr(s_hat_mid,120  );

if (  ((abs(difans))==0)==TRUE){ break};
if (  ((abs(difans))<=.00000001)==TRUE){ break};
if (  (difmaxmin<=.00000000001)==TRUE ) {
s_hat<-Rmpfr::mpfr(NA,120  );
break}
if (  (counter_inner==10000)==TRUE) {
s_hat<-Rmpfr::mpfr(NA,120  );
break}

}

ans<-Rmpfr::mpfr(s_hat,120  )
if(!is.na(ans)){break}
if ((counter_outer==5)==TRUE){
s_hat<-Rmpfr::mpfr(NA,120  );
break}
}

return(ans)
}
#' p_kprs2ECM computes the 2nd der as a function of s_hat,s_hat found by function p_find_sECM 
#' @param s_hat The value needed for the first derivative to equal log(lambdaECM),found by p_find_sECM
#' @param N1 Number of rows in matrix 1  
#' @param N2 Number of rows in matrix 2  
#' @param sigma_eigv A vector of eigenvalues of (Sigma11)*(inv-Sigma22)
p_kprs2ECM<- function(s_hat,N1,N2,sigma_eigv) {


if(is.na(s_hat)){
second_der<-NA;
} else {

parma<-10^(-8); 
parm1<-(s_hat+parma);
parm2<-(s_hat-parma);

first_der1<-(p_kprsECM(parm1,N1,N2,sigma_eigv));
first_der2<-(p_kprsECM(parm2,N1,N2,sigma_eigv));
diffr<-(first_der1)-(first_der2);
second_der<-  diffr/(2*parma);
return(second_der);
}
}
#' p_l_and_rECM computes the Lugananni and Rice approx. of power 
#' @param lambdaECM The value of a equality of covariance likelihood ratio statistic 
#' @param N1 Number of rows in matrix 1  
#' @param N2 Number of rows in matrix 2  
#' @param sigma_eigv A vector of eigenvalues of (Sigma11)*(inv-Sigma22);
p_l_and_rECM<- function(lambdaECM,N1,N2,sigma_eigv) {

s_hat<-p_find_sECM(lambdaECM,N1,N1,sigma_eigv);

if(is.na(s_hat)){
power_answers<-matrix(data=(c(NA,NA)),1,2);
} else {

lambdaECMn<-Rmpfr::mpfr(lambdaECM,precBits=120 ) 
y<-log(lambdaECMn)
k_s_hat<-p_lmglmECM(s_hat,N1,N2,sigma_eigv);
u1<-p_kprs2ECM(s_hat,N1,N2,sigma_eigv);
u<-(s_hat)*(sqrt(abs(u1)))
r1<-y*s_hat
r2<-r1-k_s_hat;
r3<-(2*(abs(r2)))^(1/2)
r<-(sign(s_hat))*r3;
mpower1<-Rmpfr::mpfr( (Rmpfr::pnorm(r,lower.tail=FALSE)) ,precBits=120    );
piv<-Rmpfr::mpfr(pi,120)
nd_r1<-Rmpfr::mpfr((1/(sqrt(2*piv))),precBits=120)
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
ans<-c(type2error);
pv_answer<-as.matrix(ans);
return(pv_answer)
}
}
#' statsECMf   Computes a lambdaECM statistics and p-values
#'
#' \code{statsECMf} takes two matricies of numeric data, the LR statistic lambdaECM, computed p values  
#' @param data_matrix1 Is a numeric matrix,
#' @param data_matrix2 Is a numeric matrix,
#' @examples
#'     suppressPackageStartupMessages(library(Rmpfr)); 
#'     library(wmpvaer);
#'     data(plastic)
#'     a1<-plastic[plastic$additive=="Low",];
#'     low_additive<-as.matrix(a1[,1:3]);
#'     a2<-plastic[!(plastic$additive=="Low"),];
#'     high_additive<-as.matrix(a2[,1:3]);
#'     statsECMf(low_additive,high_additive);
#' @keywords multivariate    &   Multivariate Techniques
#' @return One matrix with two elements  
#' @export
statsECMf<-function(data_matrix1,data_matrix2){

data_matrix_num1<-as.numeric(data_matrix1);
if( (any(is.na(data_matrix_num1)))==TRUE ){stop("no characters allowed in data_matrix1") }

data_matrix_num2<-as.numeric(data_matrix2);
if( (any(is.na(data_matrix_num2)))==TRUE ){stop("no characters allowed in data_matrix2") }

if((is.null(colnames(data_matrix1)))==TRUE){
colnames(data_matrix1)<-colnames(data_matrix,do.NULL=FALSE)
};

if((is.null(colnames(data_matrix2)))==TRUE){
colnames(data_matrix2)<-colnames(data_matrix,do.NULL=FALSE)
};

N1<-nrow(data_matrix1);
N2<-nrow(data_matrix2);
n1<-N1-1;
n2<-N2-1;
n<-n1+n2;

dm_cov1<-cov(data_matrix1);
dm_cov2<-cov(data_matrix2);

test_det_cov1<- function(covm){
if(det(covm)<=0) stop("determinant of covariance matrix1 <=0, lambdaECM can't be computed")
}
test_det_cov1(dm_cov1);

test_det_cov2<- function(covm){
if(det(covm)<=0) stop("determinant of covariance matrix2 <=0, lambdaECM can't be computed")
}
test_det_cov2(dm_cov2);

invdm_cov2<-tryCatch(solve(dm_cov2),error= function(e) stop("covariance matrix for input matrix2 not computationally invertible") )

eigmECM<- (dm_cov1 %*% invdm_cov2);
sigma_eigv<-eigen(eigmECM)$values;

a<-(det(dm_cov1))^(n1/n);
b<-(det(dm_cov2))^(n2/n);
c<-a*b;
d<-det(dm_cov1+dm_cov2);
lambdaECM<-c/d;

pw_and_p<-p_l_and_rECM(lambdaECM,N1,N2,sigma_eigv) 

if((is.na(pw_and_p[1,1]))==FALSE){
pvalue<-Rmpfr::asNumeric((pw_and_p[1,1]))
}else{
pvalue<-c("No Answer");
}

lambdaECM<-format(lambdaECM,digits=8);
pvalue<-format(pvalue,digits=8,justify=c("none"));
lanswer<-matrix(data=c(lambdaECM,pvalue),1,2);
rownames(lanswer)<-c(" ");
colnames(lanswer)<-c("lambdaECM","p value")
ans<-(noquote(lanswer))

return(ans)

}


