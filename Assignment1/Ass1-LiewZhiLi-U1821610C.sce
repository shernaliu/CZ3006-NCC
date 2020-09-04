/* define values of g from 0 to 50, step interval 0.1 */
g_val=[0:0.1:20]';

/* define values of a & p 
- Do not set a=0.01
- Do not set p=0.01 or 0.03
*/
a_val=0.02;
p_val=0.05;
q_val=1-p_val;

/* define function for pure aloha from eqn 1 */
function pure_aloha_val=pure_aloha(G)
    pure_aloha_val=G.*exp(-2.*G);
endfunction

/* define function for slotted aloha from eqn 2 */
function slotted_aloha_val=slotted_aloha(G)
    slotted_aloha_val=G.*exp(-G);
endfunction

/* define function for non-persistent CSMA from eqn 3 */
function non_persistent_csma_val=non_persistent_csma(A,G)
    numerator_np=(G.*exp(-A.*G));
    denominator_np=(G.*(1+2.*A)+exp(-A.*G));
    non_persistent_csma_val=numerator_np./denominator_np;
endfunction

/* define function for 1-persistent CSMA from eqn 10 */
function one_persistent_csma_val=one_persistent_csma(A,G)
    numerator_op=(G).*(1+G+(A.*G).*(1+G+(A.*G)./2)).*(exp(-G.*(1+2.*A)));
    denominator_op=G.*(1+2.*A)-(1-exp(-A.*G))+(1+A.*G).*(exp(-G.*(1+A))) ;
    one_persistent_csma_val=numerator_op./denominator_op;
endfunction

/* define function for p-persistent CSMA from A1 aka S' */
function p_persistent_val=p_persistent_csma(A,G,P,Q)
    numerator_pp=(1-exp(-A.*G)).*((Ps_prime_hat(A,G,P,Q).*pi_0(A,G))+(Ps_hat(A,G,P,Q).*(1-pi_0(A,G))));
    denominator_pp=(1-exp(-A.*G)).*((A.*t_bar_prime_hat(A,G).*pi_0(A,G))+(A.*t_bar_hat(A,G)).*(1-pi_0(A,G))+1+A)+(A.*pi_0(A,G));
    p_persistent_val=numerator_pp./denominator_pp;
endfunction

/* define function for pi_0 from eqn 25*/
function pi_0_val=pi_0(A,G)
    pi_0_val=exp(-G.*(1+A));
endfunction

/* define function for Ps_hat from eqn A20 using pi_0(A,G) */
function Ps_hat_val=Ps_hat(A,G,P,Q)
    g=G*A;
    eqn_a=(((pi_0(A,G)).^P)-(pi_0(A,G)))./(Q.*(1-pi_0(A,G)));
    eqn_b=((1-exp(-g.*P)).*(((pi_0(A,G)).^(1-Q.^2))-pi_0(A,G)))./((Q.*(1-pi_0(A,G)))-(Q.*exp(-2.*g.*P)).*(((pi_0(A,G)).^P)-pi_0(A,G)));
    Ps_hat_val=eqn_a-eqn_b;
endfunction

/* define function for Ps_prime_hat from eqn A20 replacing pi_0(A,G) with exp(-g) where g=aG */
function Ps_prime_hat_val=Ps_prime_hat(A,G,P,Q)
    g=G*A;
    eqn_a=(((exp(-(A.*G))).^P)-(exp(-(A.*G))))./(Q.*(1-exp(-(A.*G))));
    eqn_b=((1-exp(-g.*P)).*(((exp(-(A.*G))).^(1-Q.^2))-exp(-(A.*G))))./((Q.*(1-exp(-(A.*G))))-(Q.*exp(-2.*g.*P)).*(((exp(-(A.*G))).^P)-exp(-(A.*G))));
    Ps_prime_hat_val=eqn_a-eqn_b;
endfunction

/* define function for t-bar-hat from eqn A12 */
function t_bar_hat_val=t_bar_hat(A,G,P)
    g=G*A;
    t_bar_hat_val=((pi_0(A,G).^P)-pi_0(A,G))./(1-pi_0(A,G)-((pi_0(A,G).^P)-pi_0(A,G)).*exp(-P.*g));
endfunction

/* define function for t-bar-prime-hat from eqn A12 replacing pi-naught with exp(-g) where g=aG */
function t_bar_prime_hat_val=t_bar_prime_hat(A,G,P)
    g=G*A;
    t_bar_prime_hat_val=((exp(-(A.*G)).^P)-exp(-(A.*G)))./(1-exp(-(A.*G))-((exp(-(A.*G)).^P)-exp(-(A.*G))).*exp(-(P.*g)));
endfunction

/* --- plot2d customization options --- */
colors=[color("red") color("orange") color("yellow") color("green") color("blue")];
xtitle("Throughput Analyses for ALOHA and CSMA (a=0.02,p=0.05)");
xlabel("G (Offered Channel Traffic)","fontsize", 3,"fontname",8);
ylabel("S (Throughput)","fontsize", 3,"fontname",8);
option_a=gca();
option_a.tight_limits=["on","on"];
option_a.data_bounds=[0 0;20 1];
// change font properties
option_a.font_foreground=1; 
option_a.font_size=3;
option_a.font_style=8;
// change thickness of polyline
p=get("hdl");
p.thickness=2;

/* plot the curves on a graph */
plot2d(g_val,[pure_aloha(g_val) slotted_aloha(g_val) non_persistent_csma(a_val,g_val) one_persistent_csma(a_val,g_val) p_persistent_csma(a_val,g_val,p_val,q_val)],colors,leg="Pure ALOHA@Slotted ALOHA@Non-Persistent CSMA@1-Persistent CSMA@0.05-Persistent CSMA");

/* calculate max throughput for Pure ALOHA & Slotted ALOHA */
disp("Max. Throughput of Pure ALOHA at G=0.5 is: ", pure_aloha(0.5));
disp("Max. Throughput of Slotted ALOHA at G=1 is: ", slotted_aloha(1));

/* calculate max throughput for Non-P, 1-P, 0.05-P CSMA */
disp("Max. Throughput of Non-Persistent CSMA at G=6.5 is: ", non_persistent_csma(a_val,6.5));
disp("Max. Throughput of 1-Persistent CSMA at G=1 is: ", one_persistent_csma(a_val,1));
disp("Max. Throughput of 0.05-Persistent CSMA at G=3.2 is: ", p_persistent_csma(a_val,3.2,p_val,q_val));
