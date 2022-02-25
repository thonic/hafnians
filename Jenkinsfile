pipeline {
    agent any
    stages {
        stage('Run hafnian.py') {
            steps {
                sh '''
                    pwd
                    whoami
                    /home/thonic/anaconda3/envs/partial-derivatives/bin/python /home/thonic/git/partial-derivatives/hafnian.py
                '''
            }
        }
    }
}