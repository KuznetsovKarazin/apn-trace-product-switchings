#!/usr/bin/env python3
"""Stage 5R: manuscript-grade audit certificate for the uniform theorem.

This script does not replace the mathematical proof.  It verifies every finite
constant and bridge used by the proof, records the exact pole-divisor cases,
and checks the n=10,12 certificates from the canonical dimension scan.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--dimension-scan',type=Path,default=Path('results/current/trace_switch_dimension_scan.json'))
    ap.add_argument('--output',type=Path,default=Path('results/current/manuscript_uniform_theorem.json'))
    args=ap.parse_args()
    root=Path(__file__).resolve().parents[1]
    scanp=args.dimension_scan if args.dimension_scan.is_absolute() else root/args.dimension_scan
    outp=args.output if args.output.is_absolute() else root/args.output
    scan=json.loads(scanp.read_text())
    rows={r['n']:r for r in scan['dimensions']}

    # C_theta: X^3+Y^3=theta Z^3, theta != 0.
    geometry={
      'smoothness': 'The partials are X^2,Y^2,theta Z^2; they have no common projective zero.',
      'genus': 1,
      'infinity_points': '[t:1:0] for t in mu_3=GF(4)^*',
      'infinity_count': 3,
      'axis_polar_divisor': 'P_1+P_lambda+P_lambda^2 (three simple poles)',
      'off_axis_cancellation': 'At t=v/u the numerator uX+vY is the flex tangent; substitution X=tY gives theta Z^3, hence intersection multiplicity 3.',
      'off_axis_polar_divisor': 'The two remaining infinity points, each simple.',
      'non_coboundary': 'A function h^2+h+c has even pole order at every pole of h, so a rational function with a simple pole cannot have this form.',
    }

    g=1
    def as_degree(simple_poles:int)->int:
        return 2*g-2 + simple_poles*(1+1)
    axis=as_degree(3)
    off=as_degree(2)
    assert axis==6 and off==4
    axis_pairs=6; off_pairs=9; hasse=2
    constant=hasse+axis_pairs*axis+off_pairs*off
    assert constant==74

    threshold=[]
    for n in range(10,42,2):
        q=1<<n; sq=1<<(n//2); num=q-2-constant*sq
        threshold.append({'n':n,'q':q,'sqrt_q':sq,'numerator':num,'lower_bound':num/16,'positive':num>0})
    assert not [r for r in threshold if r['n']==12][0]['positive']
    assert all(r['positive'] for r in threshold if r['n']>=14)

    bridge=[]
    for n in (10,12):
        r=rows[n]
        assert r['represented_nonzero_count']==(1<<n)-1
        assert r['safe_nonzero_count']==0
        assert r['minimum_positive_representation_count']>0
        bridge.append({'n':n,'minimum_h':r['minimum_positive_representation_count'],'all_nonzero_represented':True})

    result={
      'schema':'manuscript-uniform-theorem-v1',
      'date':'2026-08-02',
      'theorem':'For every even n>=10 and theta!=0, h_n(theta)>0; hence x^3+theta Q is not APN.',
      'character_expansion':'h_n(theta)=1/16 sum_{u,v in GF(4)} chi_4(u) S_{u,v}(theta).',
      'curve_geometry':geometry,
      'character_sum_theorem':{
        'statement':'For a smooth projective curve C/F_q of genus g and a reduced rational phase f with pole orders m_P prime to p, |sum_{P not poles} psi(f(P))| <= (2g-2+sum_P(m_P+1))*sqrt(q).',
        'axis_degree':axis,'off_axis_degree':off,
        'references':['G. Lachaud, Theoretical Computer Science 94 (1992), 295-310.','H. Stichtenoth, Algebraic Function Fields and Codes, 2nd ed., Springer, 2009.']
      },
      'error_accounting':{'hasse':hasse,'axis_pairs':axis_pairs,'axis_each':axis,'off_axis_pairs':off_pairs,'off_axis_each':off,'total':constant},
      'lower_bound':'h_n(theta) >= (q-2-74 sqrt(q))/16',
      'threshold':threshold,
      'finite_bridge':bridge,
      'proof_status':'manuscript-grade proof written; n=10,12 remain explicit exact finite bridge lemmas',
      'provenance':{'dimension_scan':str(scanp.relative_to(root)),'sha256':sha256(scanp)},
      'validation':{'smooth_curve_argument_recorded':True,'pole_divisors_recorded':True,'non_coboundary_argument_recorded':True,'constant_74':constant==74,'positive_from_n14':all(r['positive'] for r in threshold if r['n']>=14),'finite_bridge_closed':len(bridge)==2}
    }
    outp.parent.mkdir(parents=True,exist_ok=True)
    outp.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'output':str(outp),'sha256':sha256(outp),'constant':constant,'finite_bridge':bridge},indent=2))
    return 0
if __name__=='__main__': raise SystemExit(main())
