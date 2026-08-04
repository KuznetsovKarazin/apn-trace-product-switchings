#!/usr/bin/env python3
"""Stage 5S: exact EA/CCZ partition of the eight marked projective switches.

Within each Frobenius orbit, an explicit EA witness is verified:
  G_{rho^2,eta^2}(x) = (G_{rho,eta}(x^{2^7}))^2 + L_{rho,eta}(x),
where L is an explicitly recorded F2-linear map.
The two orbits have distinct orthoderivative differential and Walsh spectra.
For quadratic APN functions this separates their EA classes and, using the
standard coincidence of EA and CCZ equivalence in the quadratic APN setting,
their CCZ classes.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
from stage5g_trace_coordinate_theory import absolute_trace,gf_inv,gf_mult,gf_pow
from stage5l_permutation_screen import build_lut
from stage5_canonical_low_rank_locus import portable_ortho_signature

N=256

def sha256(p:Path)->str: return hashlib.sha256(p.read_bytes()).hexdigest()

def is_linear(lut:list[int])->bool:
    if lut[0]!=0:return False
    return all(lut[x^y]==(lut[x]^lut[y]) for x in range(N) for y in range(N))

def linear_rank(lut:list[int])->int:
    basis=[0]*8; rank=0
    for x in [1<<i for i in range(8)]:
        v=lut[x]
        while v:
            p=v.bit_length()-1
            if basis[p]:v^=basis[p]
            else:basis[p]=v;rank+=1;break
    return rank

def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--stage5j',type=Path,default=Path('results/current/projective_displacement_theorem.json'))
    ap.add_argument('--output',type=Path,default=Path('results/current/projective_equivalence_partition.json'))
    args=ap.parse_args(); root=Path(__file__).resolve().parents[2]
    inp=args.stage5j if args.stage5j.is_absolute() else root/args.stage5j
    out=args.output if args.output.is_absolute() else root/args.output
    data=json.loads(inp.read_text())
    orbits=[]
    all_records=[]
    for oi,orbit in enumerate(data['accepted_Frobenius_orbits'],1):
        recs=[]
        for row in orbit:
            rho=int(row['rho0_hex'],16); eta=int(row['eta_equals_zeta_hex'],16)
            lut=build_lut(rho,eta); sig=portable_ortho_signature(lut)
            rec={'rho_hex':f'0x{rho:02x}','eta_hex':f'0x{eta:02x}','function_sha256':hashlib.sha256(bytes(lut)).hexdigest(),'orthoderivative_signature':sig}
            recs.append(rec); all_records.append((rho,eta,lut,sig))
        # verify every arrow in the cyclic Frobenius orbit
        witnesses=[]
        lookup={(r,e):l for r,e,l,s in all_records}
        for row in recs:
            rho=int(row['rho_hex'],16); eta=int(row['eta_hex'],16)
            trho=gf_pow(rho,2); teta=gf_pow(eta,2)
            source=build_lut(rho,eta); target=build_lut(trho,teta)
            conjugate=[gf_pow(source[gf_pow(x,128)],2) for x in range(N)]
            coeff=gf_mult(eta,gf_pow(gf_inv(rho),3))
            correction=[gf_mult(gf_pow(eta,2),absolute_trace(x)) ^ gf_mult(gf_pow(coeff,2),absolute_trace(gf_mult(gf_pow(rho,2),x))) for x in range(N)]
            assert is_linear(correction)
            assert all(target[x]==(conjugate[x]^correction[x]) for x in range(N))
            witnesses.append({'source':[f'0x{rho:02x}',f'0x{eta:02x}'],'target':[f'0x{trho:02x}',f'0x{teta:02x}'],'input_linear_map':'x -> x^(2^7)','output_linear_map':'y -> y^2','linear_correction_formula':'eta^2 Tr(x) + (eta*rho^(-3))^2 Tr(rho^2 x)','linear_correction_rank':linear_rank(correction),'identity_verified_on_256_inputs':True})
        sigs={r['orthoderivative_signature']['sha256'] for r in recs}; assert len(sigs)==1
        orbits.append({'orbit_id':oi,'members':recs,'common_signature_sha256':next(iter(sigs)),'explicit_EA_cycle':witnesses,'EA_equivalent_within_orbit':True})
    assert len(orbits)==2
    assert orbits[0]['common_signature_sha256']!=orbits[1]['common_signature_sha256']
    result={'schema':'projective-equivalence-partition-v1','date':'2026-08-02','objects':'eight marked projective APN switchings','partition':orbits,'between_orbits':{'different_orthoderivative_differential_spectra':orbits[0]['members'][0]['orthoderivative_signature']['differential_spectrum']!=orbits[1]['members'][0]['orthoderivative_signature']['differential_spectrum'],'different_orthoderivative_absolute_walsh_spectra':orbits[0]['members'][0]['orthoderivative_signature']['absolute_walsh_spectrum']!=orbits[1]['members'][0]['orthoderivative_signature']['absolute_walsh_spectrum'],'EA_inequivalent':True,'CCZ_inequivalent_for_quadratic_APN':True},'theorem':{'EA_classes':2,'CCZ_classes':2,'class_sizes_among_marked_points':[4,4],'identification':['PART2-X3-1','PART2-X3-2']},'proof_dependencies':{'within_class':'explicit EA identities verified on all field elements','between_classes':'distinct orthoderivative differential/Walsh spectra; these are EA invariants, and EA=CCZ for quadratic APN functions'},'provenance':{'stage5j':str(inp.relative_to(root)),'sha256':sha256(inp)},'validation':{'eight_points':sum(len(o['members']) for o in orbits)==8,'two_frobenius_orbits':len(orbits)==2,'all_EA_arrows_verified':all(w['identity_verified_on_256_inputs'] for o in orbits for w in o['explicit_EA_cycle']),'signatures_separate_orbits':orbits[0]['common_signature_sha256']!=orbits[1]['common_signature_sha256']}}
    out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'output':str(out),'sha256':sha256(out),'EA_classes':2,'CCZ_classes':2,'signatures':[o['common_signature_sha256'] for o in orbits]},indent=2))
    return 0
if __name__=='__main__': raise SystemExit(main())
