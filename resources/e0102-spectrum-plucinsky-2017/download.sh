wget https://wikis.mit.edu/confluence/download/attachments/60927986/NoLine_v1.3.1_coco.fits
wget https://wikis.mit.edu/confluence/download/attachments/60927986/NoLine_v1.3.1_line.fits
wget https://wikis.mit.edu/confluence/download/attachments/60927986/rgspn_mod_tbabs_tbvarabs_2apec_line_ratios_jd_v1.9.xcm
convert_xspec_script rgspn_mod_tbabs_tbvarabs_2apec_line_ratios_jd_v1.9.xcm sherpa-model-e0102.py --clobber

