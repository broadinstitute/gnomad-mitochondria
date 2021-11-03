#!/usr/bin/env python
import argparse
import os

import hail as hl
import hailtop.batch as hb

from batch.batch_utils import (
    init_arg_parser,
    init_job,
    run_batch,
    switch_gcloud_auth_to_user_account,
)


def main():  # noqa: D103
    p = init_arg_parser()
    p.add_argument(
        "--infile", required=True, help="Tab delimited file of participant and bam_path"
    )
    p.add_argument("--docker", required=True, help="Docker image to use")
    p.add_argument("--output-bucket", required=True, help="Output bucket for results")
    p.add_argument(
        "--mutserv-script",
        required=True,
        help="Path to process_mutserv.py python script",
    )
    p.add_argument("--gcloud-user", required=True, help="gcloud user account")
    p.add_argument("--gcloud-project", required=True, help="gcloud project")
    p.add_argument("--batch-job-name", help="Batch: (optional) job name")

    args = p.parse_args()

    output_bucket = args.output_bucket

    with run_batch(args) as b:

        python_script = b.read_input(args.mutserv_script)

        with hl.hadoop_open(args.infile, "r") as f:
            for line in f:
                line = line.rstrip()
                items = line.split("\t")
                participant, bam_path = items[0:2]

                output_file_path = os.path.join(output_bucket, f"{participant}.vcf.gz")
                output_split_path = os.path.join(
                    output_bucket, f"{participant}_split.vcf"
                )
                output_mutserv_path = os.path.join(
                    output_bucket, f"{participant}_split_mutserv.vcf"
                )

                file_stats = hl.hadoop_stat(bam_path)
                bam_size = int(round(file_stats["size_bytes"] / 10.0 ** 9))

                j = init_job(
                    b,
                    name=participant,
                    cpu=args.cpu,
                    memory=args.memory,
                    disk_size=bam_size * 2,
                    image=args.docker,
                )

                # Switch_gcloud_auth_to_user_account(j, args.gcloud_creds, args.gcloud_user, args.gcloud_project)
                switch_gcloud_auth_to_user_account(
                    j, args.gsa_key_file, args.gcloud_user, args.gcloud_project
                )

                # Run mutserv
                j.command(f"gsutil -m cp {bam_path} {participant}.cram")
                j.command(
                    f"java -jar mutserve-1.3.4.jar analyse-local --deletions --insertions --input {participant}.cram --output {participant}.vcf.gz --reference rCRS.fasta --level 0.01"
                )
                j.command(f"tabix -f {participant}.vcf.gz")
                # Flip alleles where so that format fields are reported in same order as the alt alleles are listed, for example, the following line in the VCF would be reformatted:
                # chrM POS . REF A,G . PASS . GT:AF:DP 2/1:0.98,0.012:2740 -> chrM POS . REF A,G . PASS . GT:AF:DP 1/2:0.012,0.98:2740
                j.command(
                    """zcat {participant}.vcf.gz | awk '{{ if ($10~/2[\/]1/) {{split($10,a,":"); split(a[2],b,","); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t1/2:"b[2]","b[1]":"a[3]}} else {{print $0}}}}' > {participant}_temp.vcf""".format(
                        **locals()
                    )
                )

                # Split multiallelics with bcftools
                j.command(
                    f"bcftools view --header-only {participant}_temp.vcf| sed 's/##FORMAT=<ID=AF,Number=1/##FORMAT=<ID=AF,Number=A/' > {participant}_pre.vcf"
                )
                j.command(
                    f"bcftools view --no-header {participant}_temp.vcf >> {participant}_pre.vcf"
                )
                j.command(
                    f"bcftools norm -m -any {participant}_pre.vcf > {participant}_split.vcf"
                )

                # Reformat output
                with open("process_mutserv.py", "rt") as f:
                    mutserv_script_code = f.read()

                j.command(
                    f"""cat > mutserv_script.py <<- EOF
{mutserv_script_code}
EOF
"""
                )
                j.command(
                    f"python3 mutserv_script.py --input-file {participant}_split.vcf --output-file {participant}_split_mutserv.txt --mt-reference rCRS.fasta"
                )

                # Copy results to bucket
                j.command(f"cp {participant}.vcf.gz {j.ofile}")
                j.command(f"cp {participant}_split.vcf {j.ofile2}")
                j.command(f"cp {participant}_split_mutserv.txt {j.ofile3}")

                b.write_output(j.ofile, output_file_path)
                b.write_output(j.ofile2, output_split_path)
                b.write_output(j.ofile3, output_mutserv_path)

        b.run(open=True)


if __name__ == "__main__":
    main()
