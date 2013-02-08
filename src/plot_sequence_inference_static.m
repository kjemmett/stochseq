function plot_sequence_inference_static(est,dna)
% function inf_obj = plot_sequence_inference_static(est,dna)
% 
% Display image object of inferred sequence against
% true DNA sequence
%
% Can be used independently or as a subroutine 
% for animation generation
% 
%
% inputs:
%   est : sequence estimate (Lx4)
%   dna : true dna sequence

figure

ax(1)=subplot(2,1,1);
inf_obj=imagesc(est');
hold on
plot(dna','y.-','LineWidth',1.5);
hold off
xlabel('Sequence Position');
ylabel('Nucleotide');
set(gca,'YTick',[1 2 3 4],'YTickLabel',{'A','C','G','T'});

ax(2)=subplot(2,1,2);
ent_obj=bar(calc_entropy(est));xlim([0 length(dna)]);ylim([0 1]);xlabel('Sequence Position');ylabel('Entropy');
set(zoom,'Motion','horizontal','Enable','on')

tot_ent = sum(calc_entropy(est));
title(['total entropy ',num2str(tot_ent)])

linkaxes([ax(1) ax(2)],'x');
