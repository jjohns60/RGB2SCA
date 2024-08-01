function [F1,precision,recall,accuracy] = F1score(A,B)
%F1score calculates the F1score from a reference binary dataset A to the
% dataset B. Inputs must have the same dimensions. The function returns the
% F1 score as well as the precision, recall, and accuracy

%prepare data
A = A(:);
B = B(:);

%calculate true positives
TP = sum(A == 1 & B == 1,"all");

%calculate false positives
FP = sum(A == 1 & B == 0,"all");

%calculate false negatives
FN = sum(A == 0 & B == 1,"all");

%calculate precision, recall, and F1
precision = (TP)/(TP + FP);
recall = (TP)/(TP + FN);
F1 = (2 * precision * recall)/(precision + recall);

%calculate the accuracy
accuracy = sum(A == B,"all")./numel(A);

end