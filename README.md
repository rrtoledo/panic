# PANIC: a Pluri Anonymous Non Interactive Credential system

This repo. includes the code of a crypto. primitive we call Group Signatures with Proof of Ownership (GSPO) bulit on on zk-proof, digital signature and cca2 encryption. This primitive can easily be leveraged into a privacy friendly identity broker.

This repo. contains the following files:
- bsps.py : A bilateral structure preserving signature ("Optimal Structure-Preserving Signature in Assymetric Bilinear Group" by Abe Groth et al.)
- cca2_cs.py : Cramer-Shoup CCA2 encryption
- cca2_eg.py : CCA2 encryption based on El Gamal and
- gspo.py : our new primitive and a running script
- gsproof.py : Groth Sahai proof ("Fine-tuning groth-sahai proofs" by Escala and Groth)
- proofs.py : the proofs needed for our GSPO primitive
- proofs_aggreg.py : the proofs needed for our GSPO primitive but aggregated
- sok.py : a Signature of Knowledge based on Fiat Shamir heuristic
- usps.py : a unilateral structure preserving signature ("Short Group Signatures via Structure-Preserving Signatures" from Libert et al.)
