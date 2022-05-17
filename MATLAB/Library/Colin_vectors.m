function [Vectors_product_2] = Colin_vectors(VTransmitter_target,VTarget_receiver)

Vectors_product_2=dot(VTransmitter_target,VTarget_receiver)/(norm(VTransmitter_target)*norm(VTarget_receiver));
angle_vectors=180-acosd(Vectors_product_2);


end