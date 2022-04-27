pipeline {
    agent any
    stages {
        stage('Run hafnian.py') {
            when { allOf {
                branch 'master';
                changelog '^RUN .*'
                }
            }
            steps {
                sh '''
                    pwd
                    m=$(git log -1 --pretty=%B | sed 's/RUN \\(.*\\)/\\1/')
                    echo $m
                    /home/thonic/anaconda3/envs/partial-derivatives/bin/python /var/lib/jenkins/workspace/hafnian_master/hafnian.py --job_name="$m"
                    git add *
                    git commit -m "save results"
                    sha=$(git rev-parse HEAD)
                    git push origin $sha:master
                '''
            }
        }
    }
}