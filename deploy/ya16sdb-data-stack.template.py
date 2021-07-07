from troposphere import Template, Ref, Tags, Join, GetAtt, Parameter, Output
from awacs.aws import Action, Allow, ArnEquals, Condition, Policy, PolicyDocument, Principal, SourceArn, Statement
from awacs import s3 as awacs_s3

import troposphere.iam as iam
import troposphere.s3  as s3


t = Template()
t.set_version()
t.set_description("Ya16sdb application stack")


t.add_parameter(
    Parameter(
        "DataBucketName",
        Description="S3 bucket name for application data",
        Type="String",
    )
)

application_policy = iam.ManagedPolicy(
    "AppYa16sdbPolicy",
    ManagedPolicyName="AppYa16sdbPolicy",
    PolicyDocument=PolicyDocument(
        Version="2012-10-17",
        Statement=[
            Statement(
                Effect=Allow,
                Action=[
                    awacs_s3.GetObject,
                    awacs_s3.PutObject,
                    awacs_s3.ListBucket,
                ],
                Resource=[
                    Join("", [ "arn:aws:s3:::", Ref("DataBucketName") ]),
                    Join("", [ "arn:aws:s3:::", Ref("DataBucketName"), "/*" ]),
                ],
            )
        ]
    )
)
t.add_resource(application_policy)


s3bucket = s3.Bucket(
    "S3Bucket",
    BucketName=Ref("DataBucketName"),
    AccessControl="Private",
    PublicAccessBlockConfiguration=s3.PublicAccessBlockConfiguration(
        BlockPublicAcls=True,
        BlockPublicPolicy=True,
        IgnorePublicAcls=True,
        RestrictPublicBuckets=True,
    ),
    BucketEncryption=s3.BucketEncryption(
        ServerSideEncryptionConfiguration=[
            s3.ServerSideEncryptionRule(
                ServerSideEncryptionByDefault=s3.ServerSideEncryptionByDefault(
                    SSEAlgorithm="AES256")
            )
        ]
    ),
    DeletionPolicy="Delete"
)
t.add_resource(s3bucket)

appuser = iam.User(
    "BucketUser",
    UserName="Ya16sdbUser",
    ManagedPolicyArns=[Ref("AppYa16sdbPolicy")]
)
t.add_resource(appuser)

print(t.to_json())