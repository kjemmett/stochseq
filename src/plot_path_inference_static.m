function plot_path_inference_static(gamma,z)
% function plot_path_inference_static(gamma,z)
%
% Return an image object specifying
%
% Input
% -----
%
% gamma : (T x L) double
%   Gamma matrix specifying p(z_t=l|X)
%   for every t.
%
% z : (T x 1) int
%   True path along the sequence. Elements will be 
%   int's from 1 to L.
%
% Output
% ------
%
% inf_z_img_obj
%   Image object containing path inference and true path 

T = size(gamma,1);
L = size(gamma,2);

ent=gamma.*log(gamma);
ent(isnan(ent)==1) = 0;
path_inf_entropy = -sum(ent,2)/log(L);

figure

ax(1)=subplot(3,2,1:4);
inf_z_img_obj = imagesc(gamma');
hold on;
plot(z,'y-');
hold off;
xlabel('Time Step');
ylabel('Sequence Position');

ax(2)=subplot(3,2,5:6);
path_inf_ent = bar(path_inf_entropy');
xlim([0 length(z)]);
ylim([0 1]);
xlabel('Time Step');
ylabel('Entropy');

linkaxes([ax(1) ax(2)],'x');
