%
%$Id: calc_equi_bw_design_space.m,v 1.4 2024/06/12 09:15:13 venkatnarayan.h Exp venkatnarayan.h $
%
% Copyright (c) 2024, Shiv Nadar University, Delhi NCR, India. All Rights
% Reserved. Permission to use, copy, modify and distribute this software for
% educational, research, and not-for-profit purposes, without fee and without a
% signed license agreement, is hereby granted, provided that this paragraph and
% the following two paragraphs appear in all copies, modifications, and
% distributions.
%
% IN NO EVENT SHALL SHIV NADAR UNIVERSITY BE LIABLE TO ANY PARTY FOR DIRECT,
% INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST
% PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE.
%
% SHIV NADAR UNIVERSITY SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT
% NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS PROVIDED "AS IS". SHIV
% NADAR UNIVERSITY HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
% ENHANCEMENTS, OR MODIFICATIONS.
%
% Revision History:
% Date           By                      Change Notes
% -------------  ----------------------  --------------------------------------
% 12th Jun 2024  V. Hariharan            Original
%

%Parameters
bw_reqd_MHz = 100;
a0_min = 100000;  %Min open loop gain spec of the basic amplifier (BA)
f = 0.005;  %Feedback factor. Ensure that a0_min.f >> 1, for approximations in analytical model to hold

m0 = 2;       %Starting value of p2/p1
mfinal = 40;  %Ending value of p2/p1
mBumpPct = 10;

n0 = 2;       %Starting value of p3/p2
nfinal = 40;  %Ending value of p3/p2
nBumpPct = 10;

disp("-I-: Starting ...\n");

%Calculated fields
T = a0_min*f;

r_m = 1 + mBumpPct/100;
count_m = 1 + log(mfinal/m0)/log(r_m);
m_vec = m0 * ((1 + mBumpPct/100) .^ (0 : count_m - 1));

p1_vec = [];
p2_vec = [];
p3_vec = [];
BW_ExactNum = [];
BW_ErrPct = [];

for m = m_vec
  r_n = 1 + nBumpPct/100;
  count_n = 1 + log(nfinal/n0)/log(r_n);
  n_vec = n0 * ((1 + nBumpPct/100) .^ (0 : count_n - 1));

  for n = n_vec
    alpha = T/m;
    beta = sqrt(1 - 4*alpha + 8*alpha^2);
    t1 = beta - 1 + 2*alpha;
    p2 = -bw_reqd_MHz * sqrt((2 + t1/n^2) / t1 / (1 + alpha^2/beta/n^2));
    p1 = p2 / m;
    p3 = p2 * n;
    exact_num_bw = calc_bandwidth(p1, p2, p3, a0_min, T);

    %Calculate err% wrt reqd BW
    bw_err_pct = 100 * abs((bw_reqd_MHz - exact_num_bw) / bw_reqd_MHz);

    p1_vec = [p1_vec, p1];
    p2_vec = [p2_vec, p2];
    p3_vec = [p3_vec, p3];
    BW_ExactNum = [BW_ExactNum, exact_num_bw];
    BW_ErrPct = [BW_ErrPct, bw_err_pct];
  end  %End loop n=p3/p2
end  %End loop m=p2/p1

max_bw_err_pct = max(BW_ErrPct);
min_bw_err_pct = min(BW_ErrPct);
avg_bw_err_pct = mean(BW_ErrPct);

csv_fname = compose("BA_PoleSpaceForBW_MHz=%.1f_T=%.2f.csv", bw_reqd_MHz, T);
writematrix(["#p1(MHz)" "p2(MHz)" "p3(MHz)" "ExactBW(MHz)"], csv_fname);
writematrix([p1_vec' p2_vec' p3_vec' BW_ExactNum'], csv_fname, ...
                                        'WriteMode', 'append');

p1_lin = linspace(min(p1_vec), max(p1_vec), 20);
p2_lin = linspace(min(p2_vec), max(p2_vec), 20);
[X, Y] = meshgrid(p1_lin, p2_lin);
F = scatteredInterpolant(p1_vec', p2_vec', p3_vec');
%F = scatteredInterpolant(p1_vec', p2_vec', p3_vec', 'nearest'); %This seems to
                                          %work better for BW=6 MHz and Tmin=80
Z = F(X, Y);

fh = figure(1);
ah = gca;

surf(ah, X, Y, Z) %interpolated

title(ah, compose("Equi-BW surface for BW=%.2f MHz and Tmin = %.2f", ...
                  bw_reqd_MHz, T));
xlabel(ah, "p1 (MHz)");
ylabel(ah, "p2 (MHz)");
zlabel(ah, "p3 (MHz)");
set(ah, 'ZDir', 'reverse');
%xlim(ah, [-0.4, 0]); %This looks better for BW=6 MHz and Tmin=80

%Save PNG image
png_fname = compose("BA_PoleSpaceForBW_MHz=%.1f_T=%.2f.png", bw_reqd_MHz, T);
exportgraphics(fh, png_fname, "Resolution", 300);
delete(fh);


%Plot the error%
fh = figure(1);
ah = gca;
epct_obj = plot(ah, BW_ErrPct, ".");

grid(ah, "on");
title(ah, compose("Abs Error%% for BW=%.2f MHz and Tmin = %.2f", ...
                  bw_reqd_MHz, T));
xlabel(ah, "Sample point index of the \{p1,p2,p3\} trio");
ylabel(ah, "Error%");
ann_txt = compose("Max=%.1f%%  Min=%.1f%%  Avg=%.1f%%", ...
                          max_bw_err_pct, min_bw_err_pct, avg_bw_err_pct);
%BEGIN: Plot specific hardcoded values
tan_h = annotation(fh, 'textbox', [0.15 0.85 0.5 0.05], 'String', ann_txt);
tan_h.FontSize = 7;
tan_h.EdgeColor = "none";
%END: Plot specific hardcoded values

%Save PNG image
png_fname = compose("AbsErr%%ForBW_MHz=%.1f_T=%.2f.png", bw_reqd_MHz, T);
exportgraphics(fh, png_fname, "Resolution", 300);
delete(fh);

disp("-I-: Done\n");


function bw = calc_bandwidth(p1, p2, p3, a0, T)
  poly1 = [-1/p1, 1];
  poly2 = [-1/p2, 1];
  poly3 = [-1/p3, 1];

  num = [a0];
  denom = [0,0,0,T] + conv(conv(poly1, poly2), poly3);
  
  A = tf(num, denom);
  bw_low_upp = bandwidth(A);
  bw = bw_low_upp(:,1);
end

